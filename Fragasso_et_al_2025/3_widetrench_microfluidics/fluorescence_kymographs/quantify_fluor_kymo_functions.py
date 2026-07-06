# -*- coding: utf-8 -*-
"""
Functions for wide-trench fluorescence kymograph analysis.

Supports any fluorescence (SYTOX, LL-37, etc.) and phase contrast channels.

Pipeline:
  1. Trench detection from phase contrast (get_trench_xycoords)
  2. Ensemble and single-trench kymograph extraction
  3. Single-trench tracking across frames
  4. Save / load outputs
  5. Pool kymographs across positions
  6. Overlay kymograph and crop+kymo plotting
  7. For LL-37-specific: quantify bottom of empty trench signal and background correction

@author: alessio fragasso
"""

import os
import glob
import pickle
import gc
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import matplotlib.patches as patches
from matplotlib.colors import LinearSegmentedColormap
from matplotlib.patches import Rectangle
from mpl_toolkits.axes_grid1 import make_axes_locatable
from scipy.signal import find_peaks
from skimage.transform import resize
import imageio.v3 as iio


'''Custom colormaps'''
black_to_green = LinearSegmentedColormap.from_list(
    "black_to_green",
    ["black", "lime"]
)

black_to_magenta = LinearSegmentedColormap.from_list(
    "black_to_magenta",
    ["black", "magenta"]
)

black_to_cyan = LinearSegmentedColormap.from_list(
    "black_to_cyan",
    ["black", "cyan"]
)


'''Processing: trench detection and crop extraction'''

def get_expanded_crop(
    img,
    y0,
    y1,
    x0,
    x1,
    pad=100,
    fill_value=np.nan,
):
    H, W = img.shape

    y0e = y0 - pad
    y1e = y1 + pad
    x0e = x0 - pad
    x1e = x1 + pad

    target_h = y1e - y0e
    target_w = x1e - x0e

    y0s = max(0, y0e)
    y1s = min(H, y1e)
    x0s = max(0, x0e)
    x1s = min(W, x1e)

    yy0 = y0s - y0e
    yy1 = yy0 + (y1s - y0s)
    xx0 = x0s - x0e
    xx1 = xx0 + (x1s - x0s)

    crop = np.full(
        (target_h, target_w),
        fill_value,
        dtype=float,
    )

    crop[yy0:yy1, xx0:xx1] = img[y0s:y1s, x0s:x1s]

    trench_rel = {
        "y0": pad,
        "y1": pad + (y1 - y0),
        "x0": pad,
        "x1": pad + (x1 - x0),
    }

    return crop, trench_rel


def bin_kymograph_y(kymo, bin_size=5):
    n_frames, n_y = kymo.shape
    n_bins = n_y // bin_size

    kymo_trim = kymo[:, :n_bins * bin_size]
    kymo_bin = kymo_trim.reshape(n_frames, n_bins, bin_size).mean(axis=2)

    return kymo_bin


def extract_expanded_trench_ensemble(
    img,
    trench_boxes,
    pad=100,
    kymo_use_trench_x_only=True,
    kymo_x_inset_px=50,
):
    """
    Extract expanded crops around each trench.
    For 2D average: use full expanded crop.
    For kymograph: average only over the inner trench x-region.
    """

    expanded_crops = []
    used_trench_boxes = []
    trench_rel_box = None

    for tr_id, (y0, y1, x0, x1) in enumerate(trench_boxes):

        result = get_expanded_crop(
            img,
            y0, y1, x0, x1,
            pad=pad,
        )

        if result is None:
            continue

        crop, trench_rel = result

        expanded_crops.append(crop)
        used_trench_boxes.append((tr_id, y0, y1, x0, x1))
        trench_rel_box = trench_rel

    if len(expanded_crops) == 0:
        return None

    expanded_crops = np.stack(expanded_crops, axis=0)
    frame_avg_2d = np.nanmean(expanded_crops, axis=0)

    if kymo_use_trench_x_only:
        x0r = trench_rel_box["x0"] + kymo_x_inset_px
        x1r = trench_rel_box["x1"] - kymo_x_inset_px

        if x1r <= x0r:
            raise ValueError(
                "kymo_x_inset_px is too large for the trench width. "
                f"Got rel_x0={trench_rel_box['x0']}, rel_x1={trench_rel_box['x1']}, "
                f"inset={kymo_x_inset_px}."
            )

        frame_y_profile = np.nanmean(frame_avg_2d[:, x0r:x1r], axis=1)

    else:
        frame_y_profile = np.nanmean(frame_avg_2d, axis=1)

    return frame_y_profile, frame_avg_2d, expanded_crops, used_trench_boxes, trench_rel_box


def get_active_x_range_from_profile(profile, padding_min_count=20):
    """
    Find x-region without padding using the dominant repeated value.
    If no strong padding value is found, return full profile range.
    """
    profile = np.asarray(profile)

    vals, counts = np.unique(profile, return_counts=True)
    pad_val = vals[np.argmax(counts)]
    pad_count = np.max(counts)

    if pad_count < padding_min_count:
        return 0, len(profile)

    valid = profile != pad_val
    valid_idx = np.where(valid)[0]

    if len(valid_idx) == 0:
        return 0, len(profile)

    x_min = int(valid_idx[0])
    x_max = int(valid_idx[-1] + 1)

    return x_min, x_max


def get_trench_xycoords(
    phase_img,
    peak_distance=60,
    peak_prominence=1000,
    trench_width_range=(70, 130),
    gap_width_range=(150, 250),
    require_neighbor_gap=True,
    edge_search_px=100,
    min_edge_drop=100,
    fixed_trench_width_px=125,
    fixed_gap_px=180,
    refined_width_tolerance_px=45,
    extrapolate_missing=True,
    trench_height_px=555,
    boundary_px=20,
    search_after=1150,
    drop_window=7,
    settle_window=40,
    show=False,
    frame=0,
    title=None,
):
    """
    Detect trench boxes as (y0, y1, x0, x1) from phase contrast image.

    1. Detect rough trench candidates between consecutive phase peaks.
    2. For each candidate, refine both edges via local minimum search.
    3. Accept only candidates with two decent edge drops and reasonable width.
    4. Fill missing trenches by extrapolating from nearby accepted trenches.
    5. Detect y1 from bottom edge; y0 = y1 - fixed trench height.
    """

    def local_min_near_peak(profile, peak_x, side, search_px, min_drop):
        peak_x = int(peak_x)
        W = len(profile)

        if side == "left":
            a = max(0, peak_x - search_px)
            b = peak_x + 1
        elif side == "right":
            a = peak_x
            b = min(W, peak_x + search_px + 1)
        else:
            raise ValueError("side must be 'left' or 'right'")

        if b <= a:
            return None

        local = profile[a:b]
        min_x = a + int(np.argmin(local))

        peak_val = profile[peak_x]
        min_val = profile[min_x]
        drop = peak_val - min_val

        if drop < min_drop:
            return None

        return int(min_x), float(drop)

    def add_unique_xcoord(xcoords, new_coord, duplicate_tol=30):
        x0_new, x1_new = new_coord
        c_new = 0.5 * (x0_new + x1_new)

        for x0, x1 in xcoords:
            c = 0.5 * (x0 + x1)
            if abs(c - c_new) < duplicate_tol:
                return xcoords

        xcoords.append(new_coord)
        return xcoords

    def extrapolate_from_good_trenches(
        good_xcoords,
        W,
        width,
        gap,
        active_x0,
        active_x1,
    ):
        if len(good_xcoords) == 0:
            return []

        pitch = width + gap
        good_xcoords = sorted(good_xcoords, key=lambda t: t[0])

        all_xcoords = []

        first_x0 = good_xcoords[0][0]
        x0 = first_x0

        left_extra = []

        while x0 - pitch >= active_x0:
            x0 -= pitch

            if x0 >= active_x0 and x0 + width <= active_x1:
                left_extra.append((int(x0), int(x0 + width)))

        for coord in reversed(left_extra):
            all_xcoords = add_unique_xcoord(all_xcoords, coord)

        for i, coord in enumerate(good_xcoords):
            x0_good, x1_good = coord

            if x0_good >= active_x0 and x1_good <= active_x1:
                all_xcoords = add_unique_xcoord(all_xcoords, coord)

            if i < len(good_xcoords) - 1:
                x0_curr = coord[0]
                x0_next = good_xcoords[i + 1][0]

                x_fill = x0_curr + pitch

                while x_fill < x0_next - 0.5 * pitch:
                    pred = (int(x_fill), int(x_fill + width))

                    if pred[0] >= active_x0 and pred[1] <= active_x1:
                        all_xcoords = add_unique_xcoord(all_xcoords, pred)

                    x_fill += pitch

        last_x0 = good_xcoords[-1][0]
        x0 = last_x0 + pitch

        while x0 + width <= active_x1:
            if x0 >= active_x0:
                all_xcoords = add_unique_xcoord(
                    all_xcoords,
                    (int(x0), int(x0 + width)),
                )

            x0 += pitch

        all_xcoords = sorted(all_xcoords, key=lambda t: t[0])
        return all_xcoords

    phase_avg_x = np.mean(phase_img, axis=0)
    H, W = phase_img.shape

    active_x0, active_x1 = get_active_x_range_from_profile(
        phase_avg_x,
        padding_min_count=20,
    )

    peaks, _ = find_peaks(
        phase_avg_x,
        distance=peak_distance,
        prominence=peak_prominence,
    )

    distances = np.diff(peaks)

    good_xcoords = []
    failed_candidates = []
    rough_candidates = []

    for i, width in enumerate(distances):

        is_trench = trench_width_range[0] <= width <= trench_width_range[1]

        left_gap = (
            i > 0
            and gap_width_range[0] <= distances[i - 1] <= gap_width_range[1]
        )

        right_gap = (
            i < len(distances) - 1
            and gap_width_range[0] <= distances[i + 1] <= gap_width_range[1]
        )

        gap_ok = (left_gap or right_gap) if require_neighbor_gap else True

        if not (is_trench and gap_ok):
            continue

        left_peak = int(peaks[i])
        right_peak = int(peaks[i + 1])

        rough_candidates.append((left_peak, right_peak))

        left_result = local_min_near_peak(
            phase_avg_x,
            left_peak,
            side="left",
            search_px=edge_search_px,
            min_drop=min_edge_drop,
        )

        right_result = local_min_near_peak(
            phase_avg_x,
            right_peak,
            side="right",
            search_px=edge_search_px,
            min_drop=min_edge_drop,
        )

        if left_result is None or right_result is None:
            failed_candidates.append((left_peak, right_peak))
            continue

        left_min_x, left_drop = left_result
        right_min_x, right_drop = right_result

        refined_width = right_min_x - left_min_x

        width_ok = (
            abs(refined_width - fixed_trench_width_px)
            <= refined_width_tolerance_px
        )

        if not width_ok:
            failed_candidates.append((left_peak, right_peak))
            continue

        center_x = 0.5 * (left_min_x + right_min_x)
        x0 = int(round(center_x - fixed_trench_width_px / 2))
        x1 = x0 + fixed_trench_width_px

        if x0 < 0 or x1 > W:
            failed_candidates.append((left_peak, right_peak))
            continue

        good_xcoords.append((x0, x1))

    if len(good_xcoords) == 0:
        print("No good individually detected trenches.")
        return [], phase_avg_x, peaks, None

    if extrapolate_missing:
        trench_xcoords = extrapolate_from_good_trenches(
            good_xcoords,
            W=W,
            width=fixed_trench_width_px,
            gap=fixed_gap_px,
            active_x0=active_x0,
            active_x1=active_x1,
        )
    else:
        trench_xcoords = sorted(good_xcoords, key=lambda t: t[0])

    if len(trench_xcoords) == 0:
        print("No final trench coordinates.")
        return [], phase_avg_x, peaks, None

    y_profiles = []

    for x0, x1 in trench_xcoords:
        crop = phase_img[:, x0:x1]
        y_profiles.append(np.mean(crop, axis=1))

    y_profiles = np.vstack(y_profiles)
    ensemble_y_profile = np.mean(y_profiles, axis=0)

    p = ensemble_y_profile
    n = len(p)

    drop_score = np.full(n, -np.inf)

    start = max(boundary_px, search_after, drop_window)
    end = n - max(boundary_px, drop_window)

    if start >= end:
        print("Invalid y-search range.")
        return [], phase_avg_x, peaks, ensemble_y_profile

    for y in range(start, end):
        before = np.mean(p[y - drop_window:y])
        after = np.mean(p[y:y + drop_window])
        drop_score[y] = before - after

    edge_idx = int(np.argmax(drop_score))

    a = edge_idx
    b = min(n, edge_idx + settle_window)

    y1 = a + int(np.argmin(p[a:b]))
    y0 = max(boundary_px, y1 - trench_height_px)

    trench_boxes = [
        (y0, y1, x0, x1)
        for x0, x1 in trench_xcoords
    ]

    if show:

        fig, (ax1, ax2) = plt.subplots(
            2, 1,
            figsize=(12, 14),
            constrained_layout=True,
        )
        if title:
            fig.suptitle(title)
        vmin, vmax = np.nanpercentile(phase_img, [5, 98])

        ax1.imshow(
            phase_img,
            cmap="gray",
            aspect="auto",
            vmin=vmin,
            vmax=vmax,
        )

        for _, _, x0, x1 in trench_boxes:
            ax1.axvline(x0, alpha=1)
            ax1.axvline(x1, alpha=1)

        ax1.axhline(y0, color="r")
        ax1.axhline(y1, color="r")
        ax1.set_title(
            f"Final trenches: {len(trench_boxes)}, frame = {frame}"
        )

        ax2.plot(phase_avg_x)
        ax2.scatter(peaks, phase_avg_x[peaks], s=25, color="r")

        for x0, x1 in rough_candidates:
            ax2.axvspan(x0, x1, alpha=0.10, color="gray")

        for x0, x1 in failed_candidates:
            ax2.axvspan(x0, x1, alpha=0.25, color="red")

        for x0, x1 in good_xcoords:
            ax2.axvspan(x0, x1, alpha=0.35, color="green")

        for x0, x1 in trench_xcoords:
            ax2.axvline(x0, color="k", alpha=0.5)
            ax2.axvline(x1, color="k", alpha=0.5)

        ax2.set_title(
            "Green = good individual detections; red = failed; black = final"
        )
        ax2.set_xlabel("x position [px]")
        ax2.axvline(active_x0, color="purple", linestyle="--")
        ax2.axvline(active_x1, color="purple", linestyle="--")

        plt.show()

        print("Peak distances:", distances)
        print("Good xcoords:", good_xcoords)
        print("Final xcoords:", trench_xcoords)
        print("y0, y1:", y0, y1)
        print("Trench height:", y1 - y0, "px")

    return trench_boxes, phase_avg_x, peaks, ensemble_y_profile


'''Save / load kymograph outputs'''

def save_channel_kymo_outputs(
    output_folder_kymo,
    pos,
    exp_name,
    AMP_name,
    rep_name,
    channel_name,
    frame_ids,
    profiles,
    avg_2d_frames,
    trench_rel_box,
    bin_size=5,
):
    if len(profiles) == 0:
        print(f"No {channel_name} profiles to save for {pos}")
        return

    if not os.path.exists(output_folder_kymo):
        os.makedirs(output_folder_kymo)

    basename = pos + "_" + exp_name + "_" + channel_name + "_kymo"

    kymo = np.vstack(profiles)
    avg_2d_frames = np.stack(avg_2d_frames, axis=0)
    avg_2d = np.mean(avg_2d_frames, axis=0)

    kymo_binned = bin_kymograph_y(kymo, bin_size=bin_size)

    np.savez_compressed(
        os.path.join(output_folder_kymo, basename + "_arrays.npz"),
        kymo=kymo,
        kymo_binned=kymo_binned,
        avg_2d=avg_2d,
        avg_2d_frames=avg_2d_frames,
        frame_ids=np.asarray(frame_ids),
    )

    kymo_df = pd.DataFrame(kymo)
    kymo_df.columns = [f"y_{i}" for i in range(kymo.shape[1])]
    kymo_df.insert(0, "frame", frame_ids)
    kymo_df.insert(0, "channel", channel_name)
    kymo_df.insert(0, "rep", rep_name)
    kymo_df.insert(0, "AMP", AMP_name)
    kymo_df.insert(0, "pos", pos)
    kymo_df.insert(0, "exp_name", exp_name)

    kymo_df.to_pickle(
        os.path.join(output_folder_kymo, basename + "_kymo_df.pkl")
    )

    metadata = {
        "pos": pos,
        "exp_name": exp_name,
        "AMP_name": AMP_name,
        "rep_name": rep_name,
        "channel_name": channel_name,
        "frame_ids": frame_ids,
        "trench_rel_box": trench_rel_box,
        "bin_size": bin_size,
        "kymo_shape": kymo.shape,
        "kymo_binned_shape": kymo_binned.shape,
        "avg2d_shape": avg_2d.shape,
    }

    with open(os.path.join(output_folder_kymo, basename + "_metadata.pkl"), "wb") as f:
        pickle.dump(metadata, f)

    print(f"Saved {channel_name} kymo outputs for {pos}: {basename}")


def load_channel_kymo_npz_files(output_folder_kymo, channel_name):
    pattern = os.path.join(
        output_folder_kymo,
        f"*_{channel_name}_kymo_arrays.npz"
    )

    npz_files = sorted(glob.glob(pattern))

    datasets = []

    for f in npz_files:
        data = np.load(f)

        metadata_path = f.replace("_arrays.npz", "_metadata.pkl")

        if os.path.exists(metadata_path):
            with open(metadata_path, "rb") as pf:
                metadata = pickle.load(pf)
        else:
            metadata = {}

        datasets.append({
            "file": f,
            "metadata": metadata,
            "frame_ids": data["frame_ids"],
            "kymo": data["kymo"],
            "kymo_binned": data["kymo_binned"],
            "avg_2d": data["avg_2d"],
            "avg_2d_frames": data["avg_2d_frames"],
        })

    print(f"Loaded {len(datasets)} {channel_name} datasets")
    return datasets


def save_trench_metadata_outputs(
    output_folder_kymo,
    pos,
    exp_name,
    trench_box_rows,
    frame_summary_rows,
):
    if not os.path.exists(output_folder_kymo):
        os.makedirs(output_folder_kymo)

    basename = pos + "_" + exp_name + "_trench_metadata"

    trench_boxes_df = pd.DataFrame(trench_box_rows)
    frame_summary_df = pd.DataFrame(frame_summary_rows)

    trench_boxes_df.to_pickle(
        os.path.join(output_folder_kymo, basename + "_trench_boxes_df.pkl")
    )

    frame_summary_df.to_pickle(
        os.path.join(output_folder_kymo, basename + "_frame_summary_df.pkl")
    )

    print(f"Saved trench metadata for {pos}")


'''Tracking trenches'''

def link_trench_records(records, max_dx=40):
    """
    Link individual trenches across frames based on x_center.
    Adds 'track_id' to each record.
    """

    if len(records) == 0:
        return records

    frames = sorted(set(r["frame"] for r in records))

    next_track_id = 0
    active_tracks = {}

    for fr in frames:
        frame_records = [r for r in records if r["frame"] == fr]
        frame_records = sorted(frame_records, key=lambda r: r["x_center"])

        used_tracks = set()

        for r in frame_records:
            x = r["x_center"]

            best_track = None
            best_dx = np.inf

            for track_id, last_x in active_tracks.items():
                if track_id in used_tracks:
                    continue

                dx = abs(x - last_x)

                if dx < best_dx:
                    best_dx = dx
                    best_track = track_id

            if best_track is not None and best_dx <= max_dx:
                r["track_id"] = best_track
                active_tracks[best_track] = x
                used_tracks.add(best_track)

            else:
                r["track_id"] = next_track_id
                active_tracks[next_track_id] = x
                used_tracks.add(next_track_id)
                next_track_id += 1

    return records


def build_single_trench_kymos(
    single_trench_records,
    trench_rel_box,
    min_frames=20,
    kymo_use_trench_x_only=True,
    kymo_x_inset_px=50,
):
    tracks = {}

    for r in single_trench_records:
        track_id = r["track_id"]

        if track_id not in tracks:
            tracks[track_id] = {
                "frames": [],
                "crops": [],
                "x_centers": [],
            }

        tracks[track_id]["frames"].append(r["frame"])
        tracks[track_id]["crops"].append(r["crop"])
        tracks[track_id]["x_centers"].append(r["x_center"])

    track_outputs = {}

    for track_id, tr in tracks.items():
        order = np.argsort(tr["frames"])

        frames = np.asarray(tr["frames"])[order]
        crops = np.stack([tr["crops"][i] for i in order], axis=0)
        x_centers = np.asarray(tr["x_centers"])[order]

        if len(frames) < min_frames:
            continue

        avg_2d = np.mean(crops, axis=0)

        if kymo_use_trench_x_only:
            x0r = trench_rel_box["x0"] + kymo_x_inset_px
            x1r = trench_rel_box["x1"] - kymo_x_inset_px

            if x1r <= x0r:
                raise ValueError(
                    f"kymo_x_inset_px={kymo_x_inset_px} is too large. "
                    f"Trench x-range is {trench_rel_box['x0']}:{trench_rel_box['x1']}."
                )

            kymo = np.mean(crops[:, :, x0r:x1r], axis=2)

        else:
            kymo = np.mean(crops, axis=2)

        track_outputs[track_id] = {
            "frames": frames,
            "crops": crops,
            "kymo": kymo,
            "avg_2d": avg_2d,
            "x_centers": x_centers,
        }

    return track_outputs


def save_single_trench_track_outputs(
    output_folder_kymo,
    pos,
    exp_name,
    AMP_name,
    rep_name,
    channel_name,
    track_outputs,
    trench_rel_box,
):
    if len(track_outputs) == 0:
        print(f"No {channel_name} single-trench tracks to save for {pos}")
        return

    if not os.path.exists(output_folder_kymo):
        os.makedirs(output_folder_kymo)

    basename = pos + "_" + exp_name + "_" + channel_name + "_single_trench_tracks"
    track_folder = os.path.join(output_folder_kymo, basename)

    if not os.path.exists(track_folder):
        os.makedirs(track_folder)

    summary_rows = []
    trajectory_rows = []

    for track_id, tr in track_outputs.items():

        frames = np.asarray(tr["frames"])
        kymo = np.asarray(tr["kymo"])
        avg_2d = np.asarray(tr["avg_2d"])
        x_centers = np.asarray(tr["x_centers"])

        np.savez_compressed(
            os.path.join(track_folder, f"track_{track_id:03d}.npz"),
            frames=frames,
            kymo=kymo,
            avg_2d=avg_2d,
            x_centers=x_centers,
        )

        summary_rows.append({
            "exp_name": exp_name,
            "AMP": AMP_name,
            "rep": rep_name,
            "pos": pos,
            "channel": channel_name,
            "track_id": track_id,
            "n_frames": len(frames),
            "first_frame": int(frames[0]),
            "last_frame": int(frames[-1]),
            "mean_x_center": float(np.mean(x_centers)),
            "std_x_center": float(np.std(x_centers)),
            "mean_signal": float(np.nanmean(kymo)),
            "max_signal": float(np.nanmax(kymo)),
        })

        for fr, xc in zip(frames, x_centers):
            trajectory_rows.append({
                "exp_name": exp_name,
                "AMP": AMP_name,
                "rep": rep_name,
                "pos": pos,
                "channel": channel_name,
                "track_id": track_id,
                "frame": int(fr),
                "x_center": float(xc),
            })

    summary_df = pd.DataFrame(summary_rows)
    trajectory_df = pd.DataFrame(trajectory_rows)

    summary_df.to_pickle(
        os.path.join(output_folder_kymo, basename + "_summary_df.pkl")
    )

    trajectory_df.to_pickle(
        os.path.join(output_folder_kymo, basename + "_trajectory_df.pkl")
    )

    with open(os.path.join(output_folder_kymo, basename + "_metadata.pkl"), "wb") as f:
        pickle.dump(
            {
                "pos": pos,
                "exp_name": exp_name,
                "AMP_name": AMP_name,
                "rep_name": rep_name,
                "channel_name": channel_name,
                "trench_rel_box": trench_rel_box,
            },
            f,
        )

    print(f"Saved {len(track_outputs)} {channel_name} single-trench tracks for {pos}")


def load_single_trench_tracks(output_folder_kymo, channel_name="sytox"):
    track_dirs = sorted(
        glob.glob(os.path.join(output_folder_kymo, f"*_{channel_name}_single_trench_tracks"))
    )

    all_tracks = []

    for track_dir in track_dirs:
        track_files = sorted(glob.glob(os.path.join(track_dir, "track_*.npz")))

        basename = os.path.basename(track_dir)

        metadata_path = os.path.join(
            output_folder_kymo,
            basename + "_metadata.pkl",
        )

        if os.path.exists(metadata_path):
            with open(metadata_path, "rb") as f:
                metadata = pickle.load(f)
        else:
            metadata = {}

        for f in track_files:
            data = np.load(f)

            track_id = int(
                os.path.basename(f)
                .replace("track_", "")
                .replace(".npz", "")
            )

            all_tracks.append({
                "track_id": track_id,
                "track_path": f,
                "track_dir": track_dir,
                "basename": basename,
                "metadata": metadata,
                "frames": data["frames"],
                "kymo": data["kymo"],
                "avg_2d": data["avg_2d"],
                "x_centers": data["x_centers"],
            })

    print(f"Loaded {len(all_tracks)} {channel_name} tracked single trenches")
    return all_tracks


'''Pooling across positions'''

def nanmean_no_warning(arr, axis=0):
    valid = np.isfinite(arr)
    count = np.sum(valid, axis=axis)

    arr0 = np.where(valid, arr, 0)
    summed = np.sum(arr0, axis=axis)

    mean = np.full_like(summed, np.nan, dtype=float)

    np.divide(
        summed,
        count,
        out=mean,
        where=count > 0,
    )

    return mean, count


def pool_channel_datasets_by_frame(
    datasets,
    ignore_zeros=False,
    return_counts=False,
):
    all_frames = sorted(
        set(np.concatenate([d["frame_ids"] for d in datasets]).astype(int))
    )

    min_kymo_y = min(d["kymo"].shape[1] for d in datasets)
    min_2d_y = min(d["avg_2d_frames"].shape[1] for d in datasets)
    min_2d_x = min(d["avg_2d_frames"].shape[2] for d in datasets)

    pooled_kymo = []
    pooled_2d_frames = []
    pooled_frame_ids = []
    n_positions_per_frame = []

    pooled_kymo_counts = []
    pooled_2d_counts = []

    for fr in all_frames:
        kymo_this_frame = []
        avg2d_this_frame = []

        for d in datasets:
            frame_ids = d["frame_ids"].astype(int)
            idx = np.where(frame_ids == fr)[0]

            if len(idx) == 0:
                continue

            idx = idx[0]

            kymo_fr = d["kymo"][idx, :min_kymo_y].astype(float)
            avg2d_fr = d["avg_2d_frames"][idx, :min_2d_y, :min_2d_x].astype(float)

            if ignore_zeros:
                kymo_fr[kymo_fr == 0] = np.nan
                avg2d_fr[avg2d_fr == 0] = np.nan

            kymo_this_frame.append(kymo_fr)
            avg2d_this_frame.append(avg2d_fr)

        if len(kymo_this_frame) == 0:
            continue

        kymo_stack = np.stack(kymo_this_frame, axis=0)
        avg2d_stack = np.stack(avg2d_this_frame, axis=0)

        kymo_mean, kymo_count = nanmean_no_warning(kymo_stack, axis=0)
        avg2d_mean, avg2d_count = nanmean_no_warning(avg2d_stack, axis=0)

        pooled_kymo.append(kymo_mean)
        pooled_2d_frames.append(avg2d_mean)

        pooled_kymo_counts.append(kymo_count)
        pooled_2d_counts.append(avg2d_count)

        pooled_frame_ids.append(fr)
        n_positions_per_frame.append(len(kymo_this_frame))

    pooled_kymo = np.vstack(pooled_kymo)
    pooled_2d_frames = np.stack(pooled_2d_frames, axis=0)

    pooled_kymo_counts = np.vstack(pooled_kymo_counts)
    pooled_2d_counts = np.stack(pooled_2d_counts, axis=0)

    if return_counts:
        return (
            np.asarray(pooled_frame_ids),
            pooled_kymo,
            pooled_2d_frames,
            np.asarray(n_positions_per_frame),
            pooled_kymo_counts,
            pooled_2d_counts,
        )

    return (
        np.asarray(pooled_frame_ids),
        pooled_kymo,
        pooled_2d_frames,
        np.asarray(n_positions_per_frame),
    )


def pool_single_tracks_by_frame(
    single_tracks,
    selected_positions=None,
    selected_track_indices=None,
    selected_pos_track_ids=None,
    ignore_zeros=True,
    min_valid_tracks=1,
    return_counts=True,
):
    """
    Pool already-computed single-trench kymos directly.

    Parameters
    ----------
    selected_positions : list[str] or None
        e.g. ["xy01", "xy03"]
    selected_track_indices : list[int] or None
        Global Python indices into single_tracks.
    selected_pos_track_ids : list[tuple[str, int]] or None
        e.g. [("xy01", 0), ("xy03", 2)] — per-position track_id.
    ignore_zeros : bool
        If True, zeros are treated as missing.
    min_valid_tracks : int
        Pixels with fewer valid single trenches are set to NaN.
    """

    chosen = []

    for i, tr in enumerate(single_tracks):
        pos = tr["metadata"].get("pos")
        track_id = int(tr["track_id"])

        keep = True

        if selected_positions is not None:
            keep = keep and (pos in selected_positions)

        if selected_track_indices is not None:
            keep = keep and (i in selected_track_indices)

        if selected_pos_track_ids is not None:
            keep = keep and ((pos, track_id) in selected_pos_track_ids)

        if keep:
            chosen.append((i, tr))

    if len(chosen) == 0:
        raise ValueError("No single tracks selected.")

    print("Pooling", len(chosen), "single-trench tracks")

    all_frames = sorted(
        set(
            np.concatenate([
                tr["frames"].astype(int)
                for _, tr in chosen
            ])
        )
    )

    min_y = min(tr["kymo"].shape[1] for _, tr in chosen)

    pooled_kymo = []
    pooled_counts = []
    pooled_frame_ids = []

    for fr in all_frames:
        kymos_this_frame = []

        for track_index, tr in chosen:
            frames = tr["frames"].astype(int)
            idx = np.where(frames == fr)[0]

            if len(idx) == 0:
                continue

            kymo_fr = tr["kymo"][idx[0], :min_y].astype(float)

            if ignore_zeros:
                kymo_fr[kymo_fr == 0] = np.nan

            kymos_this_frame.append(kymo_fr)

        if len(kymos_this_frame) == 0:
            continue

        stack = np.stack(kymos_this_frame, axis=0)

        valid = np.isfinite(stack)
        count = np.sum(valid, axis=0)

        stack0 = np.where(valid, stack, 0)
        summed = np.sum(stack0, axis=0)

        mean = np.full(stack.shape[1], np.nan, dtype=float)

        np.divide(
            summed,
            count,
            out=mean,
            where=count >= min_valid_tracks,
        )

        mean[count < min_valid_tracks] = np.nan

        pooled_kymo.append(mean)
        pooled_counts.append(count)
        pooled_frame_ids.append(fr)

    pooled_kymo = np.vstack(pooled_kymo)
    pooled_counts = np.vstack(pooled_counts)
    pooled_frame_ids = np.asarray(pooled_frame_ids)

    if return_counts:
        return pooled_frame_ids, pooled_kymo, pooled_counts, chosen

    return pooled_frame_ids, pooled_kymo


'''Overlay kymograph plotting'''

def plot_pooled_overlay_kymograph(
    pooled_frame_ids,
    pooled_kymo_signal,
    pooled_kymo_phase,
    fr_inj,
    trench_rel_box=None,
    AMP_title=None,
    signal_name="Signal",
    overlay_mode="green",       # "green" or "colormap"
    overlay_cmap="inferno",
    signal_vmin_prc=5,
    signal_vmax_prc=99,
    manual_signal_vmin=None,
    manual_signal_vmax=None,
    phase_vmin_prc=2,
    phase_vmax_prc=100,
    manual_phase_vmin=None,
    manual_phase_vmax=None,
    phase_invert=False,
    signal_alpha=0.85,
    zero_background="grey",     # "black", "white", "grey", "light_grey"
    show_trench_lines=False,
    ylims_px=None,              # e.g. (0, 700)
    save_path=None,
    px_size=0.065841,
    figsize=(9.5, 6.8),
    ytick_step_um=10,
):
    def normalize_with_limits(
        img,
        vmin_prc=2,
        vmax_prc=99,
        manual_vmin=None,
        manual_vmax=None,
        ignore_zeros=True,
    ):
        img = img.astype(float)

        if ignore_zeros:
            vals = img[np.isfinite(img) & (img != 0)]
        else:
            vals = img[np.isfinite(img)]

        if len(vals) == 0:
            return np.zeros_like(img, dtype=float), 0.0, 1.0

        if manual_vmin is not None and manual_vmax is not None:
            vmin = float(manual_vmin)
            vmax = float(manual_vmax)
        else:
            vmin, vmax = np.nanpercentile(vals, [vmin_prc, vmax_prc])

        norm = (img - vmin) / (vmax - vmin + 1e-9)
        norm = np.clip(norm, 0, 1)

        return norm, float(vmin), float(vmax)

    signal_kymo = pooled_kymo_signal.astype(float)
    phase_kymo = pooled_kymo_phase.astype(float)

    if signal_kymo.shape != phase_kymo.shape:
        raise ValueError(
            f"Signal and phase pooled kymos must have same shape, got "
            f"{signal_kymo.shape} and {phase_kymo.shape}"
        )

    time_min = pooled_frame_ids - fr_inj

    phase_norm, phase_vmin, phase_vmax = normalize_with_limits(
        phase_kymo,
        vmin_prc=phase_vmin_prc,
        vmax_prc=phase_vmax_prc,
        manual_vmin=manual_phase_vmin,
        manual_vmax=manual_phase_vmax,
        ignore_zeros=True,
    )

    signal_norm, signal_vmin, signal_vmax = normalize_with_limits(
        signal_kymo,
        vmin_prc=signal_vmin_prc,
        vmax_prc=signal_vmax_prc,
        manual_vmin=manual_signal_vmin,
        manual_vmax=manual_signal_vmax,
        ignore_zeros=True,
    )

    if phase_invert:
        phase_norm = 1 - phase_norm

    phase_invalid = (~np.isfinite(phase_kymo)) | (phase_kymo == 0)
    signal_invalid = ~np.isfinite(signal_kymo)

    phase_disp = phase_norm.T
    signal_disp = signal_norm.T

    phase_invalid_disp = phase_invalid.T
    signal_invalid_disp = signal_invalid.T

    if zero_background == "black":
        bg_val = 0.0
    elif zero_background == "white":
        bg_val = 1.0
    elif zero_background in ["grey", "gray"]:
        bg_val = 0.5
    elif zero_background in ["light_grey", "light_gray"]:
        bg_val = 0.85
    else:
        raise ValueError(
            "zero_background must be 'black', 'white', 'grey', or 'light_grey'."
        )

    phase_disp[phase_invalid_disp] = bg_val
    signal_disp[signal_invalid_disp] = 0
    signal_disp[phase_invalid_disp] = 0

    rgb = np.zeros((*phase_disp.shape, 3), dtype=float)
    rgb[..., 0] = phase_disp
    rgb[..., 1] = phase_disp
    rgb[..., 2] = phase_disp

    rgb[~np.isfinite(rgb)] = bg_val
    signal_disp[~np.isfinite(signal_disp)] = 0

    alpha = signal_alpha * signal_disp

    if overlay_mode == "green":
        rgb[..., 0] = rgb[..., 0] * (1 - alpha)
        rgb[..., 1] = rgb[..., 1] * (1 - alpha) + alpha
        rgb[..., 2] = rgb[..., 2] * (1 - alpha)

        signal_cmap_for_cbar = LinearSegmentedColormap.from_list(
            "signal_green",
            ["black", "lime"]
        )

    elif overlay_mode == "colormap":
        cmap_obj = plt.get_cmap(overlay_cmap)
        signal_rgba = cmap_obj(signal_disp)
        signal_rgb = signal_rgba[..., :3]
        rgb = rgb * (1 - alpha[..., None]) + signal_rgb * alpha[..., None]
        signal_cmap_for_cbar = cmap_obj

    else:
        raise ValueError("overlay_mode must be 'green' or 'colormap'.")

    fig, ax = plt.subplots(figsize=figsize)
    fig.subplots_adjust(right=0.86)

    fig.patch.set_facecolor("white")
    ax.set_facecolor((0.85, 0.85, 0.85))

    y0_um = 0.0
    y1_um = signal_kymo.shape[1] * px_size

    ax.imshow(
        rgb,
        aspect="auto",
        origin="upper",
        extent=[
            time_min[0],
            time_min[-1],
            y1_um,
            y0_um,
        ],
    )

    if show_trench_lines and trench_rel_box is not None:
        line_color = "white" if zero_background != "white" else "black"

        ax.axhline(
            trench_rel_box["y0"],
            color=line_color,
            linestyle="--",
            linewidth=1.5,
        )
        ax.axhline(
            trench_rel_box["y1"],
            color=line_color,
            linestyle="--",
            linewidth=1.5,
        )

    y_max_um = signal_kymo.shape[1] * px_size

    if ylims_px is not None:
        y_start_px, y_end_px = ylims_px

        y_start_um = y_start_px * px_size
        y_end_um = y_end_px * px_size

        ax.set_ylim(
            y_max_um - y_start_um,
            y_max_um - y_end_um,
        )

        visible_span_um = y_end_um - y_start_um

        ytick_labels_um = np.arange(
            0,
            visible_span_um + ytick_step_um,
            ytick_step_um,
        )

        ytick_positions = y_max_um - (y_start_um + ytick_labels_um)

    else:
        ytick_labels_um = np.arange(
            0,
            y_max_um + ytick_step_um,
            ytick_step_um,
        )
        ytick_positions = y_max_um - ytick_labels_um

    ylim_low, ylim_high = ax.get_ylim()
    ymin_vis = min(ylim_low, ylim_high)
    ymax_vis = max(ylim_low, ylim_high)

    keep = (ytick_positions >= ymin_vis) & (ytick_positions <= ymax_vis)

    ax.set_yticks(ytick_positions[keep])
    ax.set_yticklabels([f"{y:.0f}" for y in ytick_labels_um[keep]])
    ax.tick_params(axis="both", labelsize=14)

    ax.axvline(x=0, linestyle=":", c="black")

    ax.set_xlabel("Time since start of AMP addition (min)", fontsize=14)
    ax.set_ylabel("Position along trench (μm)", fontsize=14)

    title = f"Pooled {signal_name}/phase overlay kymograph"
    if AMP_title is not None:
        title = AMP_title + "\n" + title
    ax.set_title(title, fontsize=13)

    cax_signal = fig.add_axes([0.88, 0.15, 0.02, 0.30])
    sm_signal = plt.cm.ScalarMappable(
        cmap=signal_cmap_for_cbar,
        norm=plt.Normalize(vmin=signal_vmin, vmax=signal_vmax),
    )
    sm_signal.set_array([])
    cbar_signal = fig.colorbar(sm_signal, cax=cax_signal)
    cbar_signal.set_label(f"{signal_name} intensity", fontsize=11)

    if save_path is not None:
        fig.savefig(save_path, format="pdf", dpi=300, bbox_inches="tight")

    plt.show()


def get_dataset_by_pos(datasets, pos):
    for d in datasets:
        if d["metadata"].get("pos") == pos:
            return d
    raise ValueError(f"No dataset found for pos={pos}")


def plot_overlay_kymo_for_one_position(
    output_folder_kymo,
    pos,
    fr_inj,
    channel_name="sytox",
    AMP_title=None,
    signal_name="Signal",
    save_path=None,
    overlay_mode="colormap",
    overlay_cmap="inferno",
    manual_signal_vmin=None,
    manual_signal_vmax=None,
    signal_vmin_prc=10,
    signal_vmax_prc=90,
    phase_vmin_prc=15,
    phase_vmax_prc=95,
    phase_invert=False,
    signal_alpha=0.85,
    zero_background="light_grey",
    show_trench_lines=False,
    ylims_px=(250, 1200),
    ytick_step_um=10,
):
    signal_datasets = load_channel_kymo_npz_files(
        output_folder_kymo,
        channel_name=channel_name,
    )

    phase_datasets = load_channel_kymo_npz_files(
        output_folder_kymo,
        channel_name="phase",
    )

    ds_signal = get_dataset_by_pos(signal_datasets, pos)
    ds_phase = get_dataset_by_pos(phase_datasets, pos)

    frame_ids_signal = ds_signal["frame_ids"].astype(int)
    frame_ids_phase = ds_phase["frame_ids"].astype(int)

    if not np.array_equal(frame_ids_signal, frame_ids_phase):
        raise ValueError(
            f"Frame IDs do not match for pos={pos}: "
            f"{frame_ids_signal.shape} vs {frame_ids_phase.shape}"
        )

    trench_rel_box = ds_signal["metadata"].get("trench_rel_box", None)

    plot_pooled_overlay_kymograph(
        pooled_frame_ids=frame_ids_signal,
        pooled_kymo_signal=ds_signal["kymo"],
        pooled_kymo_phase=ds_phase["kymo"],
        fr_inj=fr_inj,
        trench_rel_box=trench_rel_box,
        AMP_title=(f"{AMP_title}\n{pos}" if AMP_title is not None else pos),
        signal_name=signal_name,
        overlay_mode=overlay_mode,
        overlay_cmap=overlay_cmap,
        manual_signal_vmin=manual_signal_vmin,
        manual_signal_vmax=manual_signal_vmax,
        signal_vmin_prc=signal_vmin_prc,
        signal_vmax_prc=signal_vmax_prc,
        phase_vmin_prc=phase_vmin_prc,
        phase_vmax_prc=phase_vmax_prc,
        phase_invert=phase_invert,
        signal_alpha=signal_alpha,
        zero_background=zero_background,
        show_trench_lines=show_trench_lines,
        ylims_px=ylims_px,
        ytick_step_um=ytick_step_um,
        save_path=save_path,
    )


'''Single-trench crop + kymograph plotting (helpers)'''

def list_files_like_processing(path, use_sorted=False):
    files = [f for f in os.listdir(path) if ".DS_Store" not in f]
    return sorted(files) if use_sorted else files


def load_trench_boxes_df_for_track(output_folder_kymo, track):
    meta = track["metadata"]

    pos = meta["pos"]
    exp_name = meta["exp_name"]

    pattern = os.path.join(
        output_folder_kymo,
        f"{pos}_{exp_name}_trench_metadata_trench_boxes_df.pkl"
    )

    matches = glob.glob(pattern)

    if len(matches) == 0:
        raise FileNotFoundError(f"No trench metadata file found for pattern:\n{pattern}")

    trench_boxes_df = pd.read_pickle(matches[0])

    trench_boxes_df["x_center"] = 0.5 * (
        trench_boxes_df["x0"] + trench_boxes_df["x1"]
    )

    return trench_boxes_df


def get_expanded_crop_asym(
    img,
    y0,
    y1,
    x0,
    x1,
    pad_x=100,
    pad_top=250,
    pad_bottom=100,
):
    H, W = img.shape

    y0e = max(0, y0 - pad_top)
    y1e = min(H, y1 + pad_bottom)
    x0e = max(0, x0 - pad_x)
    x1e = min(W, x1 + pad_x)

    crop = img[y0e:y1e, x0e:x1e]

    trench_rel = {
        "y0": y0 - y0e,
        "y1": y1 - y0e,
        "x0": x0 - x0e,
        "x1": x1 - x0e,
    }

    return crop, trench_rel


def resize_image_to_match(img, target_shape, order=1):
    if img.shape == target_shape:
        return img.astype(float)

    img_resized = resize(
        img,
        target_shape,
        order=order,
        preserve_range=True,
        anti_aliasing=True,
    )

    return img_resized.astype(float)


def extract_tracked_trench_crops_for_frames(
    single_tracks_sytox,
    track_index,
    frames_to_extract,
    output_folder_kymo,
    master_folder_path,
    pad=100,
    max_frame_delta=2,
    fluor_channel="fluor2",
    use_sorted=False,
    pad_x=100,
    pad_top=250,
    pad_bottom=100,
    signal_bkg_sub=False,
):
    """
    Re-read raw images and extract selected tracked-trench crops.
    Phase crop is raw (non-background-subtracted).
    Fluorescence crop is optionally background-subtracted.
    The crop is matched to the tracked trench by x_center.
    """

    track = single_tracks_sytox[track_index]
    meta = track["metadata"]

    AMP_name = meta["AMP_name"]
    rep_name = meta["rep_name"]
    exp_name = meta["exp_name"]
    pos = meta["pos"]

    experiment_path = os.path.join(master_folder_path, AMP_name, rep_name)
    folder_path = os.path.join(experiment_path, pos)

    mask_path = os.path.join(folder_path, "masks")
    phase_path = os.path.join(folder_path, "phase")
    fluor_path = os.path.join(folder_path, fluor_channel)

    masks_list = list_files_like_processing(mask_path, use_sorted=use_sorted)
    phase_list = list_files_like_processing(phase_path, use_sorted=use_sorted)
    fluor_list = list_files_like_processing(fluor_path, use_sorted=use_sorted)

    trench_boxes_df = load_trench_boxes_df_for_track(
        output_folder_kymo,
        track,
    )

    track_frames = track["frames"].astype(int)
    track_xcenters = track["x_centers"]

    crop_records = []

    for target_frame in frames_to_extract:

        idx_track = int(np.argmin(np.abs(track_frames - target_frame)))
        matched_frame = int(track_frames[idx_track])

        if abs(matched_frame - target_frame) > max_frame_delta:
            print(
                f"Skipping target frame {target_frame}: "
                f"closest tracked frame is {matched_frame}"
            )
            continue

        track_x_center = float(track_xcenters[idx_track])

        df_fr = trench_boxes_df[trench_boxes_df["frame"] == matched_frame].copy()

        if len(df_fr) == 0:
            print(f"No trench boxes found for frame {matched_frame}")
            continue

        df_fr["dx_to_track"] = np.abs(df_fr["x_center"] - track_x_center)
        row = df_fr.sort_values("dx_to_track").iloc[0]

        y0 = int(row["y0"])
        y1 = int(row["y1"])
        x0 = int(row["x0"])
        x1 = int(row["x1"])

        mask_temp = iio.imread(os.path.join(mask_path, masks_list[matched_frame]))
        phase_temp = iio.imread(os.path.join(phase_path, phase_list[matched_frame]))
        fluor_temp = iio.imread(os.path.join(fluor_path, fluor_list[matched_frame]))

        if signal_bkg_sub:
            bs = int(fluor_temp.shape[0] / 16)

            fluor_signal = local_bkg_sub_cp(
                fluor_temp,
                mask_temp,
                pos,
                exp_name,
                frame=str(matched_frame),
                box_size=bs,
                dilations=15,
                sigma_=60,
                show=False,
            )[0]

        else:
            fluor_signal = fluor_temp.copy()

        phase_result = get_expanded_crop_asym(
            phase_temp,
            y0, y1, x0, x1,
            pad_x=pad_x,
            pad_top=pad_top,
            pad_bottom=pad_bottom,
        )

        sytox_result = get_expanded_crop_asym(
            fluor_signal,
            y0, y1, x0, x1,
            pad_x=pad_x,
            pad_top=pad_top,
            pad_bottom=pad_bottom,
        )

        if phase_result is None or sytox_result is None:
            print(f"Skipping frame {matched_frame}: crop hits image boundary")
            continue

        phase_crop_raw, phase_trench_rel = phase_result
        sytox_crop_bkg_sub, sytox_trench_rel = sytox_result

        if phase_crop_raw.shape != sytox_crop_bkg_sub.shape:
            phase_crop_raw = resize_image_to_match(
                phase_crop_raw,
                sytox_crop_bkg_sub.shape,
                order=1,
            )

        trench_rel = sytox_trench_rel

        crop_records.append({
            "target_frame": target_frame,
            "frame": matched_frame,
            "time_min": matched_frame,
            "track_index": track_index,
            "track_id": track["track_id"],
            "pos": pos,
            "AMP": AMP_name,
            "rep": rep_name,
            "exp_name": exp_name,
            "x_center_track": track_x_center,
            "x_center_box": float(row["x_center"]),
            "dx_to_track": float(row["dx_to_track"]),
            "y0": y0,
            "y1": y1,
            "x0": x0,
            "x1": x1,
            "pad": pad,
            "trench_rel_box": trench_rel,
            "phase_crop_raw": phase_crop_raw,
            "sytox_crop_bkg_sub": sytox_crop_bkg_sub,
            "pad_x": pad_x,
            "pad_top": pad_top,
            "pad_bottom": pad_bottom,
        })

    return crop_records


def _normalize_with_limits(
    img,
    vmin=None,
    vmax=None,
    vmin_prc=2,
    vmax_prc=99,
    ignore_zeros=True,
):
    img = img.astype(float)

    if ignore_zeros:
        vals = img[np.isfinite(img) & (img != 0)]
    else:
        vals = img[np.isfinite(img)]

    if len(vals) == 0:
        return np.zeros_like(img, dtype=float), 0.0, 1.0

    if vmin is None or vmax is None:
        vmin, vmax = np.nanpercentile(vals, [vmin_prc, vmax_prc])

    out = (img - vmin) / (vmax - vmin + 1e-9)
    out = np.clip(out, 0, 1)

    return out, float(vmin), float(vmax)


def _make_phase_signal_overlay(
    phase_img,
    signal_img,
    signal_vmin,
    signal_vmax,
    phase_vmin=None,
    phase_vmax=None,
    phase_vmin_prc=2,
    phase_vmax_prc=100,
    phase_invert=False,
    signal_alpha=0.85,
    zero_background="light_grey",
    overlay_mode="green",       # "green" or "colormap"
    overlay_cmap="inferno",
):
    phase_norm, _, _ = _normalize_with_limits(
        phase_img,
        vmin=phase_vmin,
        vmax=phase_vmax,
        vmin_prc=phase_vmin_prc,
        vmax_prc=phase_vmax_prc,
        ignore_zeros=True,
    )

    if phase_invert:
        phase_norm = 1 - phase_norm

    signal_norm, _, _ = _normalize_with_limits(
        signal_img,
        vmin=signal_vmin,
        vmax=signal_vmax,
        ignore_zeros=True,
    )

    if zero_background == "black":
        bg_val = 0.0
    elif zero_background == "white":
        bg_val = 1.0
    elif zero_background in ["grey", "gray"]:
        bg_val = 0.5
    elif zero_background in ["light_grey", "light_gray"]:
        bg_val = 0.85
    else:
        raise ValueError(
            "zero_background must be 'black', 'white', 'grey', or 'light_grey'."
        )

    zero_mask = phase_img == 0
    phase_norm[zero_mask] = bg_val

    rgb = np.zeros((*phase_norm.shape, 3), dtype=float)
    rgb[..., 0] = phase_norm
    rgb[..., 1] = phase_norm
    rgb[..., 2] = phase_norm

    alpha = signal_alpha * signal_norm

    if overlay_mode == "green":
        rgb[..., 0] = rgb[..., 0] * (1 - alpha)
        rgb[..., 1] = rgb[..., 1] * (1 - alpha) + alpha
        rgb[..., 2] = rgb[..., 2] * (1 - alpha)

    elif overlay_mode == "colormap":
        cmap_obj = plt.get_cmap(overlay_cmap)
        signal_rgb = cmap_obj(signal_norm)[..., :3]
        rgb = rgb * (1 - alpha[..., None]) + signal_rgb * alpha[..., None]

    else:
        raise ValueError("overlay_mode must be 'green' or 'colormap'.")

    return rgb


def _apply_relative_y_axis(
    ax,
    n_y_px,
    px_size=0.065841,
    ylims_px=None,
    ytick_step_um=10,
):
    """
    Label y-axis in um with 0 at ylims_px[0], increasing downward.
    If ylims_px=(50, 800), 50 px from the top is labeled 0.
    """
    y_max_um = n_y_px * px_size

    if ylims_px is not None:
        y_start_px, y_end_px = ylims_px

        y_start_um = y_start_px * px_size
        y_end_um = y_end_px * px_size

        ax.set_ylim(
            y_max_um - y_start_um,
            y_max_um - y_end_um,
        )

        visible_span_um = y_end_um - y_start_um

        ytick_labels_um = np.arange(
            0,
            visible_span_um + ytick_step_um,
            ytick_step_um,
        )

        ytick_positions = y_max_um - (y_start_um + ytick_labels_um)

    else:
        ytick_labels_um = np.arange(
            0,
            y_max_um + ytick_step_um,
            ytick_step_um,
        )

        ytick_positions = y_max_um - ytick_labels_um

    ylim_low, ylim_high = ax.get_ylim()
    ymin_vis = min(ylim_low, ylim_high)
    ymax_vis = max(ylim_low, ylim_high)

    keep = (ytick_positions >= ymin_vis) & (ytick_positions <= ymax_vis)

    ax.set_yticks(ytick_positions[keep])
    ax.set_yticklabels([f"{y:.0f}" for y in ytick_labels_um[keep]])


def align_kymo_to_crop_y(
    kymo,
    crop_h,
    kymo_trench_rel,
    crop_trench_rel,
    anchor="y1",
    fill_value=0,
):
    """
    Pad/shift kymo to match the height of the crop, aligning by trench position.
    anchor='y1' aligns trench bottoms; anchor='y0' aligns trench tops.
    """
    kymo = np.asarray(kymo)
    n_frames, kymo_h = kymo.shape

    out = np.full(
        (n_frames, crop_h),
        fill_value,
        dtype=kymo.dtype,
    )

    shift = int(round(crop_trench_rel[anchor] - kymo_trench_rel[anchor]))

    src0 = max(0, -shift)
    src1 = min(kymo_h, crop_h - shift)

    dst0 = max(0, shift)
    dst1 = dst0 + (src1 - src0)

    if src1 > src0:
        out[:, dst0:dst1] = kymo[:, src0:src1]

    return out, shift


def find_track_indices_by_pos(single_tracks, pos):
    return [
        i for i, tr in enumerate(single_tracks)
        if tr["metadata"].get("pos") == pos
    ]


def find_track_by_pos_and_track_id(single_tracks, pos, track_id):
    matches = [
        i for i, tr in enumerate(single_tracks)
        if (
            tr["metadata"].get("pos") == pos
            and int(tr["track_id"]) == int(track_id)
        )
    ]

    if len(matches) == 0:
        raise ValueError(f"No track found for pos={pos}, track_id={track_id}")

    if len(matches) > 1:
        print("Warning: multiple matches found:", matches)

    return matches[0]


def safe_filename(s):
    s = str(s)
    bad_chars = ['\\', '/', ':', '*', '?', '"', '<', '>', '|', ' ', ',', '\n']
    for ch in bad_chars:
        s = s.replace(ch, "_")
    return s


def get_track_match_key(track):
    meta = track["metadata"]

    return (
        meta.get("exp_name", ""),
        meta.get("pos", ""),
        int(track["track_id"]),
    )


def make_phase_track_lookup(single_tracks_phase):
    phase_lookup = {}

    for i, tr in enumerate(single_tracks_phase):
        key = get_track_match_key(tr)
        phase_lookup[key] = i

    return phase_lookup


def plot_selected_crop_and_track_kymo(
    crop_records,
    single_tracks_sytox,
    single_tracks_phase,
    crop_record_index=0,
    track_index=0,
    fr_inj=0,
    AMP_title=None,
    signal_name="SYTOX",
    px_size=0.065841,

    # signal contrast
    signal_vmin_prc=1,
    signal_vmax_prc=95,
    manual_signal_vmin=None,
    manual_signal_vmax=None,

    # phase contrast (independent for crop and kymo)
    crop_phase_vmin_prc=5,
    crop_phase_vmax_prc=100,
    kymo_phase_vmin_prc=1,
    kymo_phase_vmax_prc=100,
    manual_crop_phase_vmin=None,
    manual_crop_phase_vmax=None,
    manual_kymo_phase_vmin=None,
    manual_kymo_phase_vmax=None,

    phase_invert=False,
    signal_alpha=0.85,
    zero_background="light_grey",
    overlay_mode="green",
    overlay_cmap="inferno",
    scale_bar_um=5,

    # display
    kymo_tlims_min=None,        # e.g. (0, 250)
    ylims_px=None,              # e.g. (50, 800), 0 label starts at 50 px
    ytick_step_um=10,

    figsize=(12, 5.8),
    save_path=None,
    show=True,
):
    rec = crop_records[crop_record_index]

    phase_crop = rec["phase_crop_raw"]
    sytox_crop = rec["sytox_crop_bkg_sub"]
    frame_show = int(rec["frame"])
    time_show_min = frame_show - fr_inj

    tr_s = single_tracks_sytox[track_index]
    tr_p = single_tracks_phase[track_index]

    frames_s = tr_s["frames"].astype(int)
    frames_p = tr_p["frames"].astype(int)

    common_frames = np.intersect1d(frames_s, frames_p)

    if len(common_frames) == 0:
        raise ValueError("No common frames between signal and phase tracks.")

    idx_s = np.array([np.where(frames_s == fr)[0][0] for fr in common_frames])
    idx_p = np.array([np.where(frames_p == fr)[0][0] for fr in common_frames])

    sytox_kymo = tr_s["kymo"][idx_s]
    phase_kymo = tr_p["kymo"][idx_p]
    time_min = common_frames - fr_inj

    crop_h = phase_crop.shape[0]

    crop_trench_rel = rec["trench_rel_box"]
    kymo_trench_rel = tr_s["metadata"].get("trench_rel_box", None)

    if kymo_trench_rel is None:
        raise ValueError("Missing trench_rel_box in kymo metadata.")

    sytox_kymo, y_shift_px = align_kymo_to_crop_y(
        sytox_kymo,
        crop_h=crop_h,
        kymo_trench_rel=kymo_trench_rel,
        crop_trench_rel=crop_trench_rel,
        anchor="y1",
        fill_value=0,
    )

    phase_kymo, _ = align_kymo_to_crop_y(
        phase_kymo,
        crop_h=crop_h,
        kymo_trench_rel=kymo_trench_rel,
        crop_trench_rel=crop_trench_rel,
        anchor="y1",
        fill_value=0,
    )

    print("Applied kymo y-shift:", y_shift_px, "px")

    if manual_signal_vmin is not None and manual_signal_vmax is not None:
        signal_vmin = float(manual_signal_vmin)
        signal_vmax = float(manual_signal_vmax)

    else:
        crop_vals = sytox_crop[np.isfinite(sytox_crop) & (sytox_crop != 0)]
        kymo_vals = sytox_kymo[np.isfinite(sytox_kymo) & (sytox_kymo != 0)]

        vals = np.concatenate([crop_vals.ravel(), kymo_vals.ravel()])

        signal_vmin, signal_vmax = np.nanpercentile(
            vals,
            [signal_vmin_prc, signal_vmax_prc],
        )

        signal_vmin = float(signal_vmin)
        signal_vmax = float(signal_vmax)

    crop_rgb = _make_phase_signal_overlay(
        phase_crop,
        sytox_crop,
        signal_vmin=signal_vmin,
        signal_vmax=signal_vmax,
        phase_vmin=manual_crop_phase_vmin,
        phase_vmax=manual_crop_phase_vmax,
        phase_vmin_prc=crop_phase_vmin_prc,
        phase_vmax_prc=crop_phase_vmax_prc,
        phase_invert=phase_invert,
        signal_alpha=signal_alpha,
        zero_background=zero_background,
        overlay_mode=overlay_mode,
        overlay_cmap=overlay_cmap,
    )

    kymo_rgb = _make_phase_signal_overlay(
        phase_kymo.T,
        sytox_kymo.T,
        signal_vmin=signal_vmin,
        signal_vmax=signal_vmax,
        phase_vmin=manual_kymo_phase_vmin,
        phase_vmax=manual_kymo_phase_vmax,
        phase_vmin_prc=kymo_phase_vmin_prc,
        phase_vmax_prc=kymo_phase_vmax_prc,
        phase_invert=phase_invert,
        signal_alpha=signal_alpha,
        zero_background=zero_background,
        overlay_mode=overlay_mode,
        overlay_cmap=overlay_cmap,
    )

    fig, (ax0, ax1) = plt.subplots(
        1,
        2,
        figsize=figsize,
        gridspec_kw={
            "width_ratios": [0.55, 1.8],
            "wspace": 0.25,
        },
    )

    fig.patch.set_facecolor("white")

    h, w = phase_crop.shape
    crop_y_max_um = h * px_size

    ax0.imshow(
        crop_rgb,
        origin="upper",
        aspect="equal",
        extent=[
            0,
            w * px_size,
            crop_y_max_um,
            0,
        ],
    )

    _apply_relative_y_axis(
        ax0,
        n_y_px=h,
        px_size=px_size,
        ylims_px=ylims_px,
        ytick_step_um=ytick_step_um,
    )

    ax0.set_xlabel("x position (μm)", fontsize=12)
    ax0.set_ylabel("Position along trench (μm)", fontsize=12)
    ax0.set_title(f"t = {time_show_min} min", fontsize=13)
    ax0.tick_params(axis="both", labelsize=11)
    ax0.set_xlim(0, w * px_size)

    x_margin_um = 1.0
    y_margin_um = 1.5

    xlim = ax0.get_xlim()
    ylim = ax0.get_ylim()

    x_right = max(xlim) - x_margin_um
    x_left = x_right - scale_bar_um
    y_bottom = max(ylim) - y_margin_um

    ax0.plot(
        [x_left, x_right],
        [y_bottom, y_bottom],
        color="white",
        linewidth=3,
        solid_capstyle="butt",
    )

    ax0.text(
        0.5 * (x_left + x_right),
        y_bottom - 1.2,
        f"{scale_bar_um:g} µm",
        color="white",
        ha="center",
        va="bottom",
        fontsize=11,
    )

    kymo_y_max_um = crop_h * px_size

    ax1.imshow(
        kymo_rgb,
        origin="upper",
        aspect="auto",
        extent=[
            time_min[0],
            time_min[-1],
            kymo_y_max_um,
            0,
        ],
    )

    if kymo_tlims_min is not None:
        ax1.set_xlim(kymo_tlims_min[0], kymo_tlims_min[1])

    _apply_relative_y_axis(
        ax1,
        n_y_px=sytox_kymo.shape[1],
        px_size=px_size,
        ylims_px=ylims_px,
        ytick_step_um=ytick_step_um,
    )

    ax1.axvline(time_show_min, linestyle="--", color="white", linewidth=1.2)

    ax1.set_xlabel("Time since start of AMP addition (min)", fontsize=12)
    ax1.set_ylabel("Position along trench (μm)", fontsize=12)
    ax1.tick_params(axis="both", labelsize=11)

    title = f"{signal_name}/phase overlay"
    if AMP_title is not None:
        title = AMP_title

    fig.suptitle(title, fontsize=15)

    if overlay_mode == "green":
        signal_cmap = LinearSegmentedColormap.from_list(
            "signal_green",
            ["black", "lime"],
        )
    elif overlay_mode == "colormap":
        signal_cmap = plt.get_cmap(overlay_cmap)
    else:
        raise ValueError("overlay_mode must be 'green' or 'colormap'.")

    sm_signal = plt.cm.ScalarMappable(
        cmap=signal_cmap,
        norm=plt.Normalize(vmin=signal_vmin, vmax=signal_vmax),
    )

    sm_signal.set_array([])

    divider = make_axes_locatable(ax1)
    cax = divider.append_axes("right", size="4%", pad=0.08)

    cbar = fig.colorbar(sm_signal, cax=cax)
    cbar.set_label(f"{signal_name} intensity", fontsize=12)
    cbar.ax.tick_params(labelsize=11)

    if save_path is not None:
        fig.savefig(save_path, format="pdf", dpi=300, bbox_inches="tight")

    if show:
        plt.show()
    else:
        plt.close(fig)


def batch_plot_all_track_crop_kymos(
    single_tracks_sytox,
    single_tracks_phase,
    frames_to_extract,
    output_folder_kymo,
    master_folder_path,
    output_folder_fig,
    fr_inj,
    AMP_title=None,
    signal_name="SYTOX",

    # selection filters
    selected_positions=None,        # e.g. ["xy01", "xy04"]
    selected_track_ids=None,        # e.g. [0, 3, 7]
    selected_pos_track_pairs=None,  # e.g. [("xy01", 3), ("xy04", 7)]

    # crop extraction
    fluor_channel="fluor1",
    pad=100,
    pad_x=100,
    pad_top=250,
    pad_bottom=100,
    max_frame_delta=2,

    # signal contrast
    signal_vmin_prc=5,
    signal_vmax_prc=95,
    manual_signal_vmin=None,
    manual_signal_vmax=None,

    # phase contrast
    crop_phase_vmin_prc=5,
    crop_phase_vmax_prc=100,
    kymo_phase_vmin_prc=1,
    kymo_phase_vmax_prc=100,

    phase_invert=False,
    signal_alpha=0.85,
    zero_background="light_grey",

    overlay_mode="green",
    overlay_cmap="inferno",
    scale_bar_um=5,

    kymo_tlims_min=(0, 250),
    ylims_px=(50, 800),
    ytick_step_um=10,

    px_size=0.065841,
    use_sorted=False,
    overwrite=False,
    show=False,
):
    if not os.path.exists(output_folder_fig):
        os.makedirs(output_folder_fig)

    if selected_positions is not None:
        selected_positions = set(selected_positions)

    if selected_track_ids is not None:
        selected_track_ids = set(selected_track_ids)

    if selected_pos_track_pairs is not None:
        selected_pos_track_pairs = set(
            [(str(pos), int(tid)) for pos, tid in selected_pos_track_pairs]
        )

    phase_lookup = make_phase_track_lookup(single_tracks_phase)

    saved_rows = []

    for sytox_track_index, tr_s in enumerate(single_tracks_sytox):

        key = get_track_match_key(tr_s)

        if key not in phase_lookup:
            print(f"Skipping unmatched signal track: {key}")
            continue

        phase_track_index = phase_lookup[key]

        meta = tr_s["metadata"]
        pos = meta.get("pos", "posNA")
        exp_name = meta.get("exp_name", "expNA")
        track_id = int(tr_s["track_id"])

        if selected_positions is not None and pos not in selected_positions:
            continue

        if selected_track_ids is not None and track_id not in selected_track_ids:
            continue

        if selected_pos_track_pairs is not None and (pos, track_id) not in selected_pos_track_pairs:
            continue

        print(f"Processing {pos}, track {track_id}")

        crop_records = extract_tracked_trench_crops_for_frames(
            single_tracks_sytox=single_tracks_sytox,
            track_index=sytox_track_index,
            frames_to_extract=frames_to_extract,
            output_folder_kymo=output_folder_kymo,
            master_folder_path=master_folder_path,
            fluor_channel=fluor_channel,
            pad=pad,
            pad_x=pad_x,
            pad_top=pad_top,
            pad_bottom=pad_bottom,
            max_frame_delta=max_frame_delta,
            use_sorted=use_sorted,
        )

        if len(crop_records) == 0:
            print(f"No crop records extracted for {pos}, track {track_id}")
            continue

        for crop_i, rec in enumerate(crop_records):

            frame_used = int(rec["frame"])
            time_since_amp = frame_used - fr_inj

            filename = (
                f"{safe_filename(exp_name)}_"
                f"{safe_filename(pos)}_"
                f"track_{track_id:03d}_"
                f"frame_{frame_used:04d}_"
                f"t_{time_since_amp:+04d}min_"
                f"crop_kymo_overlay.pdf"
            )

            save_path = os.path.join(output_folder_fig, filename)

            if os.path.exists(save_path) and not overwrite:
                print(f"Already exists, skipping: {filename}")
                continue

            plot_selected_crop_and_track_kymo(
                crop_records=crop_records,
                single_tracks_sytox=[tr_s],
                single_tracks_phase=[single_tracks_phase[phase_track_index]],
                crop_record_index=crop_i,
                track_index=0,
                fr_inj=fr_inj,
                AMP_title=AMP_title,
                signal_name=signal_name,
                px_size=px_size,

                signal_vmin_prc=signal_vmin_prc,
                signal_vmax_prc=signal_vmax_prc,
                manual_signal_vmin=manual_signal_vmin,
                manual_signal_vmax=manual_signal_vmax,

                crop_phase_vmin_prc=crop_phase_vmin_prc,
                crop_phase_vmax_prc=crop_phase_vmax_prc,
                kymo_phase_vmin_prc=kymo_phase_vmin_prc,
                kymo_phase_vmax_prc=kymo_phase_vmax_prc,

                phase_invert=phase_invert,
                signal_alpha=signal_alpha,
                zero_background=zero_background,

                overlay_mode=overlay_mode,
                overlay_cmap=overlay_cmap,
                scale_bar_um=scale_bar_um,

                kymo_tlims_min=kymo_tlims_min,
                ylims_px=ylims_px,
                ytick_step_um=ytick_step_um,

                save_path=save_path,
                show=show,
            )

            saved_rows.append({
                "exp_name": exp_name,
                "pos": pos,
                "track_id": track_id,
                "frame": frame_used,
                "time_since_amp_min": time_since_amp,
                "save_path": save_path,
                "x_center_track": rec.get("x_center_track", np.nan),
                "x_center_box": rec.get("x_center_box", np.nan),
                "dx_to_track": rec.get("dx_to_track", np.nan),
            })

    saved_df = pd.DataFrame(saved_rows)

    saved_df.to_pickle(
        os.path.join(output_folder_fig, "saved_crop_kymo_overlay_index_df.pkl")
    )

    print(f"Saved {len(saved_df)} figures to:")
    print(output_folder_fig)

    return saved_df


'''LL-37 specific: bottom-trench signal and background correction'''

def plot_bottom_roi_on_crop(
    crop_records,
    crop_record_index=0,
    AMP_title=None,
    signal_name="LL-37",
    px_size=0.065841,

    manual_signal_vmin=95,
    manual_signal_vmax=140,
    crop_phase_vmin_prc=5,
    crop_phase_vmax_prc=95,
    phase_invert=False,
    signal_alpha=1,
    zero_background="black",

    overlay_mode="colormap",
    overlay_cmap="inferno",

    bottom_band_px=40,
    bottom_offset_px=0,
    kymo_x_inset_px=50,

    ylims_px=(45, 1000),
    ytick_step_um=10,
    scale_bar_um=5,

    roi_color="cyan",
    roi_linewidth=2,
    save_path=None,
    show=True,
):
    """Show a single crop with the bottom-trench ROI rectangle overlaid."""

    rec = crop_records[crop_record_index]

    phase_crop = rec["phase_crop_raw"]
    signal_crop = rec["sytox_crop_bkg_sub"]
    trench_rel_box = rec["trench_rel_box"]
    frame_show = int(rec["frame"])

    crop_rgb = _make_phase_signal_overlay(
        phase_crop,
        signal_crop,
        signal_vmin=manual_signal_vmin,
        signal_vmax=manual_signal_vmax,
        phase_vmin_prc=crop_phase_vmin_prc,
        phase_vmax_prc=crop_phase_vmax_prc,
        phase_invert=phase_invert,
        signal_alpha=signal_alpha,
        zero_background=zero_background,
        overlay_mode=overlay_mode,
        overlay_cmap=overlay_cmap,
    )

    x0_roi = max(0, int(trench_rel_box["x0"] + kymo_x_inset_px))
    x1_roi = min(signal_crop.shape[1], int(trench_rel_box["x1"] - kymo_x_inset_px))
    y1_roi = min(signal_crop.shape[0], int(trench_rel_box["y1"] - bottom_offset_px))
    y0_roi = max(0, int(y1_roi - bottom_band_px))

    fig, ax = plt.subplots(figsize=(4.5, 7))
    fig.patch.set_facecolor("white")

    h, w = phase_crop.shape
    ax.imshow(
        crop_rgb,
        origin="upper",
        aspect="equal",
        extent=[0, w * px_size, h * px_size, 0],
    )
    ax.set_xlim(0, w * px_size)

    _apply_relative_y_axis(
        ax,
        n_y_px=h,
        px_size=px_size,
        ylims_px=ylims_px,
        ytick_step_um=ytick_step_um,
    )

    rect = Rectangle(
        (x0_roi * px_size, y0_roi * px_size),
        (x1_roi - x0_roi) * px_size,
        (y1_roi - y0_roi) * px_size,
        fill=False,
        edgecolor=roi_color,
        linewidth=roi_linewidth,
    )
    ax.add_patch(rect)

    x_right = max(ax.get_xlim()) - 1.0
    x_left = x_right - scale_bar_um
    y_bottom = max(ax.get_ylim()) - 1.5

    ax.plot([x_left, x_right], [y_bottom, y_bottom],
            color="white", linewidth=3, solid_capstyle="butt")
    ax.text(0.5 * (x_left + x_right), y_bottom - 1.2,
            f"{scale_bar_um:g} µm", color="white",
            ha="center", va="bottom", fontsize=11)

    ax.set_xlabel("x position (μm)", fontsize=12)
    ax.set_ylabel("Position along trench (μm)", fontsize=12)

    title = f"{signal_name} / phase crop with bottom ROI\nframe {frame_show}"
    if AMP_title is not None:
        title = AMP_title + "\n" + title
    ax.set_title(title, fontsize=13)
    ax.tick_params(axis="both", labelsize=11)

    plt.tight_layout()

    if save_path is not None:
        os.makedirs(os.path.dirname(save_path), exist_ok=True)
        fig.savefig(save_path, format="pdf", dpi=300, bbox_inches="tight")

    if show:
        plt.show()
    else:
        plt.close(fig)

    return {"x0_roi": x0_roi, "x1_roi": x1_roi, "y0_roi": y0_roi, "y1_roi": y1_roi}


def plot_bottom_roi_crops_and_corrected_traces(
    crop_records_example,
    crop_record_indices,
    bottom_signal_corrected_df,
    fr_inj=0,
    AMP_title=None,
    signal_name="LL-37",
    px_size=0.065841,

    bottom_band_px=40,
    bottom_offset_px=0,
    kymo_x_inset_px=50,

    low_DR=True,
    crop_phase_vmin_prc=5,
    crop_phase_vmax_prc=95,
    phase_invert=False,

    ylims_px=(45, 1000),
    ytick_step_um=10,
    scale_bar_um=5,

    roi_color="cyan",
    roi_linewidth=2,
    trace_linewidth=2,
    figsize=(14, 8),

    save_path=None,
    show=True,
    xlims=None,
):
    """
    Row of crop panels (with bottom ROI box) + corrected signal trace below.
    Top: one panel per index in crop_record_indices.
    Bottom: full-width background-corrected signal vs. time.
    """

    if low_DR:
        manual_signal_vmin = 95
        manual_signal_vmax = 140
        zero_background = "black"
        signal_alpha = 1
        DR = "low_DR"
    else:
        manual_signal_vmin = 80
        manual_signal_vmax = 350
        zero_background = "light_grey"
        signal_alpha = 0.9
        DR = "high_DR"

    df_trace = bottom_signal_corrected_df.copy()

    if "signal_corrected_raw_like" not in df_trace.columns:
        raise ValueError(
            "bottom_signal_corrected_df must contain 'signal_corrected_raw_like'. "
            "Run plot_raw_like_corrected_ratio() first."
        )

    n_crops = len(crop_record_indices)

    fig = plt.figure(figsize=figsize)
    gs = fig.add_gridspec(
        2, n_crops,
        height_ratios=[1.0, 0.9],
        hspace=0.35,
        wspace=0.25,
    )

    crop_axes = [fig.add_subplot(gs[0, i]) for i in range(n_crops)]
    ax_trace = fig.add_subplot(gs[1, :])
    fig.patch.set_facecolor("white")

    shown_times = []

    for ax, crop_i in zip(crop_axes, crop_record_indices):
        rec = crop_records_example[crop_i]

        phase_crop = rec["phase_crop_raw"]
        signal_crop = rec["sytox_crop_bkg_sub"]
        trench_rel_box = rec["trench_rel_box"]
        frame_show = int(rec["frame"])
        time_show_min = frame_show - fr_inj
        shown_times.append(time_show_min)

        crop_rgb = _make_phase_signal_overlay(
            phase_crop,
            signal_crop,
            signal_vmin=manual_signal_vmin,
            signal_vmax=manual_signal_vmax,
            phase_vmin_prc=crop_phase_vmin_prc,
            phase_vmax_prc=crop_phase_vmax_prc,
            phase_invert=phase_invert,
            signal_alpha=signal_alpha,
            zero_background=zero_background,
            overlay_mode="colormap",
            overlay_cmap="inferno",
        )

        h, w = phase_crop.shape
        ax.imshow(
            crop_rgb,
            origin="upper",
            aspect="equal",
            extent=[0, w * px_size, h * px_size, 0],
        )
        ax.set_xlim(0, w * px_size)

        _apply_relative_y_axis(
            ax,
            n_y_px=h,
            px_size=px_size,
            ylims_px=ylims_px,
            ytick_step_um=ytick_step_um,
        )

        x0_roi = max(0, int(trench_rel_box["x0"] + kymo_x_inset_px))
        x1_roi = min(signal_crop.shape[1], int(trench_rel_box["x1"] - kymo_x_inset_px))
        y1_roi = min(signal_crop.shape[0], int(trench_rel_box["y1"] - bottom_offset_px))
        y0_roi = max(0, int(y1_roi - bottom_band_px))

        rect = Rectangle(
            (x0_roi * px_size, y0_roi * px_size),
            (x1_roi - x0_roi) * px_size,
            (y1_roi - y0_roi) * px_size,
            fill=False,
            edgecolor=roi_color,
            linewidth=roi_linewidth,
        )
        ax.add_patch(rect)

        x_right = max(ax.get_xlim()) - 1.0
        x_left = x_right - scale_bar_um
        y_bottom = max(ax.get_ylim()) - 1.5

        ax.plot([x_left, x_right], [y_bottom, y_bottom],
                color="white", linewidth=3, solid_capstyle="butt")
        ax.text(0.5 * (x_left + x_right), y_bottom - 1.2,
                f"{scale_bar_um:g} µm", color="white",
                ha="center", va="bottom", fontsize=10)

        ax.set_title(f"t = {time_show_min} min", fontsize=12)
        ax.set_xlabel("x (μm)", fontsize=11)

        if ax is crop_axes[0]:
            ax.set_ylabel("Position along trench (μm)", fontsize=11)
        else:
            ax.set_ylabel("")

        ax.tick_params(axis="both", labelsize=10)

    sm_signal = plt.cm.ScalarMappable(
        cmap=plt.get_cmap("inferno"),
        norm=plt.Normalize(vmin=manual_signal_vmin, vmax=manual_signal_vmax),
    )
    sm_signal.set_array([])
    divider = make_axes_locatable(crop_axes[-1])
    cax = divider.append_axes("right", size="5%", pad=0.08)
    cbar = fig.colorbar(sm_signal, cax=cax)
    cbar.set_label(f"{signal_name} intensity", fontsize=11)
    cbar.ax.tick_params(labelsize=10)

    for label in df_trace["label"].unique():
        df_tr = df_trace[df_trace["label"] == label].sort_values("time_min")
        ax_trace.plot(
            df_tr["time_min"],
            df_tr["signal_corrected_raw_like"],
            linewidth=trace_linewidth,
            label=label,
        )

    ax_trace.axvline(0, linestyle=":", color="black", linewidth=1.2)

    for tt in shown_times:
        ax_trace.axvline(tt, linestyle="--", color="gray", linewidth=1.0, alpha=0.8)

    ax_trace.set_xlabel("Time since start of AMP addition (min)", fontsize=12)
    ax_trace.set_ylabel(f"Background-corrected {signal_name} intensity", fontsize=12)
    ax_trace.tick_params(axis="both", labelsize=11)
    ax_trace.legend(frameon=False, fontsize=10)

    if xlims is not None:
        ax_trace.set_xlim(xlims)

    title = f"Background-corrected {signal_name} ROI at bottom of trench ({DR})"
    if AMP_title is not None:
        title = AMP_title + "\n" + title
    fig.suptitle(title, fontsize=14)

    plt.tight_layout()

    if save_path is not None:
        os.makedirs(os.path.dirname(save_path), exist_ok=True)
        fig.savefig(save_path, format="pdf", dpi=300, bbox_inches="tight")

    if show:
        plt.show()
    else:
        plt.close(fig)

    return df_trace


def bottom_trench_signal_from_track(
    tr,
    fr_inj=0,
    bottom_band_px=40,
    bottom_offset_px=0,
):
    """
    Compute mean signal at the bottom of the trench over time from the saved kymo.
    ROI: y1 - bottom_offset_px - bottom_band_px : y1 - bottom_offset_px
    """

    frames = tr["frames"].astype(int)
    time_min = frames - fr_inj
    kymo = tr["kymo"].astype(float)

    trench_rel_box = tr["metadata"].get("trench_rel_box", None)

    if trench_rel_box is None:
        raise ValueError("Track metadata does not contain trench_rel_box.")

    y1 = int(trench_rel_box["y1"])

    y_roi_1 = y1 - bottom_offset_px
    y_roi_0 = y_roi_1 - bottom_band_px

    y_roi_0 = max(0, y_roi_0)
    y_roi_1 = min(kymo.shape[1], y_roi_1)

    if y_roi_1 <= y_roi_0:
        raise ValueError(
            f"Invalid bottom ROI: y_roi_0={y_roi_0}, y_roi_1={y_roi_1}"
        )

    signal = np.nanmean(kymo[:, y_roi_0:y_roi_1], axis=1)

    return pd.DataFrame({
        "frame": frames,
        "time_min": time_min,
        "signal": signal,
        "track_id": int(tr["track_id"]),
        "pos": tr["metadata"].get("pos"),
        "y_roi_0": y_roi_0,
        "y_roi_1": y_roi_1,
        "bottom_y1": y1,
        "bottom_band_px": bottom_band_px,
        "bottom_offset_px": bottom_offset_px,
    })


def plot_bottom_trench_signals_for_pairs(
    single_tracks_signal,
    pos_track_pairs,
    fr_inj=0,
    bottom_band_px=40,
    bottom_offset_px=0,
    AMP_title=None,
    signal_name="LL-37",
    save_path=None,
):
    dfs = []

    for pos, track_id in pos_track_pairs:
        track_index = find_track_by_pos_and_track_id(
            single_tracks_signal,
            pos=pos,
            track_id=track_id,
        )

        tr = single_tracks_signal[track_index]

        df = bottom_trench_signal_from_track(
            tr,
            fr_inj=fr_inj,
            bottom_band_px=bottom_band_px,
            bottom_offset_px=bottom_offset_px,
        )

        df["track_index"] = track_index
        df["label"] = f"{pos}, track {track_id}"
        dfs.append(df)

    df_all = pd.concat(dfs, ignore_index=True)

    fig, ax = plt.subplots(figsize=(7, 5))
    fig.patch.set_facecolor("white")

    for label in df_all["label"].unique():
        df_tr = df_all[df_all["label"] == label]

        ax.plot(
            df_tr["time_min"],
            df_tr["signal"],
            linewidth=2,
            label=label,
        )

    ax.axvline(0, linestyle=":", color="black", linewidth=1.2)

    ax.set_xlabel("Time since start of AMP addition (min)", fontsize=13)
    ax.set_ylabel(f"Mean {signal_name} signal at trench bottom", fontsize=13)

    title = f"{signal_name} signal at bottom of trench"
    if AMP_title is not None:
        title = AMP_title + "\n" + title

    ax.set_title(title, fontsize=14)
    ax.legend(frameon=False)
    ax.tick_params(axis="both", labelsize=12)

    plt.tight_layout()

    if save_path is not None:
        os.makedirs(os.path.dirname(save_path), exist_ok=True)
        fig.savefig(save_path, format="pdf", dpi=300, bbox_inches="tight")

    plt.show()

    return df_all


def bottom_trench_signal_ratio_from_track(
    tr,
    output_folder_kymo,
    master_folder_path,
    fr_inj=0,
    fluor_channel="fluor2",
    frames_to_use=None,
    use_sorted=False,
    bottom_band_px=40,
    bottom_offset_px=0,
    kymo_x_inset_px=50,
    bg_offset_from_center_px=150,
    bg_width_px=None,
    signal_bkg_sub=False,
):
    """
    Compute signal at the bottom of the trench divided by a local background
    measured to the right of the same trench (corrects for LED flickering).
    Returns dataframe with signal_mean, background_mean, signal_bg_ratio.
    """

    meta = tr["metadata"]

    AMP_name = meta["AMP_name"]
    rep_name = meta["rep_name"]
    exp_name = meta["exp_name"]
    pos = meta["pos"]

    experiment_path = os.path.join(master_folder_path, AMP_name, rep_name)
    folder_path = os.path.join(experiment_path, pos)

    mask_path = os.path.join(folder_path, "masks")
    fluor_path = os.path.join(folder_path, fluor_channel)

    masks_list = list_files_like_processing(mask_path, use_sorted=use_sorted)
    fluor_list = list_files_like_processing(fluor_path, use_sorted=use_sorted)

    trench_boxes_df = load_trench_boxes_df_for_track(output_folder_kymo, tr)

    frames_track = tr["frames"].astype(int)
    xcenters_track = tr["x_centers"].astype(float)

    if frames_to_use is None:
        frames_to_use = frames_track
    else:
        frames_to_use = np.asarray(frames_to_use).astype(int)

    rows = []

    for target_frame in frames_to_use:

        idx_track = int(np.argmin(np.abs(frames_track - target_frame)))
        matched_frame = int(frames_track[idx_track])
        track_x_center = float(xcenters_track[idx_track])

        df_fr = trench_boxes_df[trench_boxes_df["frame"] == matched_frame].copy()

        if len(df_fr) == 0:
            print(f"No trench boxes for frame {matched_frame}")
            continue

        df_fr["dx_to_track"] = np.abs(df_fr["x_center"] - track_x_center)
        row = df_fr.sort_values("dx_to_track").iloc[0]

        y0 = int(row["y0"])
        y1 = int(row["y1"])
        x0 = int(row["x0"])
        x1 = int(row["x1"])

        fluor_temp = iio.imread(
            os.path.join(fluor_path, fluor_list[matched_frame])
        )

        if signal_bkg_sub:
            mask_temp = iio.imread(
                os.path.join(mask_path, masks_list[matched_frame])
            )

            bs = int(fluor_temp.shape[0] / 16)

            fluor_signal = local_bkg_sub_cp(
                fluor_temp,
                mask_temp,
                pos,
                exp_name,
                frame=str(matched_frame),
                box_size=bs,
                dilations=15,
                sigma_=60,
                show=False,
            )[0]
        else:
            fluor_signal = fluor_temp.astype(float)

        H, W = fluor_signal.shape

        sig_y1 = y1 - bottom_offset_px
        sig_y0 = sig_y1 - bottom_band_px

        sig_y0 = max(0, sig_y0)
        sig_y1 = min(H, sig_y1)

        sig_x0 = x0 + kymo_x_inset_px
        sig_x1 = x1 - kymo_x_inset_px

        if sig_x1 <= sig_x0:
            xc = int(round(0.5 * (x0 + x1)))
            sig_x0 = max(0, xc - 5)
            sig_x1 = min(W, xc + 5)

        sig_x0 = max(0, sig_x0)
        sig_x1 = min(W, sig_x1)

        if bg_width_px is None:
            bg_width_px_this = sig_x1 - sig_x0
        else:
            bg_width_px_this = bg_width_px

        trench_center_x = 0.5 * (x0 + x1)
        bg_center_x = trench_center_x + bg_offset_from_center_px

        bg_x0 = int(round(bg_center_x - bg_width_px_this / 2))
        bg_x1 = int(round(bg_center_x + bg_width_px_this / 2))

        bg_x0 = max(0, bg_x0)
        bg_x1 = min(W, bg_x1)

        bg_y0 = sig_y0
        bg_y1 = sig_y1

        if sig_y1 <= sig_y0 or sig_x1 <= sig_x0:
            print(f"Invalid signal ROI at frame {matched_frame}")
            continue

        if bg_y1 <= bg_y0 or bg_x1 <= bg_x0:
            print(f"Invalid background ROI at frame {matched_frame}")
            continue

        signal_mean = float(np.nanmean(fluor_signal[sig_y0:sig_y1, sig_x0:sig_x1]))
        background_mean = float(np.nanmean(fluor_signal[bg_y0:bg_y1, bg_x0:bg_x1]))

        signal_bg_ratio = signal_mean / background_mean if background_mean != 0 else np.nan

        rows.append({
            "frame": matched_frame,
            "time_min": matched_frame - fr_inj,
            "pos": pos,
            "track_id": int(tr["track_id"]),
            "signal_mean": signal_mean,
            "background_mean": background_mean,
            "signal_bg_ratio": signal_bg_ratio,
            "sig_y0": sig_y0,
            "sig_y1": sig_y1,
            "sig_x0": sig_x0,
            "sig_x1": sig_x1,
            "bg_y0": bg_y0,
            "bg_y1": bg_y1,
            "bg_x0": bg_x0,
            "bg_x1": bg_x1,
            "bottom_band_px": bottom_band_px,
            "bottom_offset_px": bottom_offset_px,
            "kymo_x_inset_px": kymo_x_inset_px,
            "bg_offset_from_center_px": bg_offset_from_center_px,
            "bg_width_px": bg_width_px_this,
        })

    return pd.DataFrame(rows)


def plot_bottom_trench_signal_ratios_for_track_indices(
    single_tracks_signal,
    track_indices,
    output_folder_kymo,
    master_folder_path,
    fr_inj=0,
    fluor_channel="fluor2",
    bottom_band_px=40,
    bottom_offset_px=0,
    kymo_x_inset_px=50,
    bg_offset_from_center_px=150,
    bg_width_px=None,
    signal_bkg_sub=False,
    AMP_title=None,
    signal_name="LL-37",
    save_path=None,
):
    dfs = []

    for track_index in track_indices:
        tr = single_tracks_signal[track_index]

        df = bottom_trench_signal_ratio_from_track(
            tr,
            output_folder_kymo=output_folder_kymo,
            master_folder_path=master_folder_path,
            fr_inj=fr_inj,
            fluor_channel=fluor_channel,
            bottom_band_px=bottom_band_px,
            bottom_offset_px=bottom_offset_px,
            kymo_x_inset_px=kymo_x_inset_px,
            bg_offset_from_center_px=bg_offset_from_center_px,
            bg_width_px=bg_width_px,
            signal_bkg_sub=signal_bkg_sub,
        )

        pos = tr["metadata"].get("pos")
        track_id = int(tr["track_id"])

        df["track_index"] = track_index
        df["label"] = f"{pos}, track {track_id} (idx {track_index})"

        dfs.append(df)

    df_all = pd.concat(dfs, ignore_index=True)

    fig, axes = plt.subplots(
        3, 1,
        figsize=(7, 8),
        sharex=True,
        constrained_layout=True,
    )

    fig.patch.set_facecolor("white")

    for label in df_all["label"].unique():
        d = df_all[df_all["label"] == label].sort_values("time_min")

        axes[0].plot(d["time_min"], d["signal_mean"], linewidth=2, label=label)
        axes[1].plot(d["time_min"], d["background_mean"], linewidth=2, label=label)
        axes[2].plot(d["time_min"], d["signal_bg_ratio"], linewidth=2, label=label)

    for ax in axes:
        ax.axvline(0, linestyle=":", color="black", linewidth=1.2)
        ax.tick_params(axis="both", labelsize=11)

    axes[0].set_ylabel("Raw signal", fontsize=12)
    axes[1].set_ylabel("Background", fontsize=12)
    axes[2].set_ylabel("Signal / background", fontsize=12)
    axes[2].set_xlabel("Time since start of AMP addition (min)", fontsize=13)

    axes[0].legend(frameon=False, fontsize=9)

    title = f"{signal_name} bottom ROI and local background"
    if AMP_title is not None:
        title = AMP_title + "\n" + title

    fig.suptitle(title, fontsize=14)

    if save_path is not None:
        os.makedirs(os.path.dirname(save_path), exist_ok=True)
        fig.savefig(save_path, format="pdf", dpi=300, bbox_inches="tight")

    plt.show()

    return df_all


def plot_raw_like_corrected_ratio(
    bottom_signal_ratio_df,
    AMP_title=None,
    signal_name="LL-37",
    ref_time_window=(-np.inf, 0),
    save_path=None,
):
    """
    Corrected signal anchored to pre-injection raw level:
        corrected = (signal / background) / mean_pre(signal / background) * mean_pre(signal)
    """

    df = bottom_signal_ratio_df.copy()

    if "signal_bg_ratio" not in df.columns:
        df["signal_bg_ratio"] = df["signal_mean"] / df["background_mean"]

    df["signal_corrected_raw_like"] = np.nan

    t0, t1 = ref_time_window

    for label in df["label"].unique():
        m = df["label"] == label
        ref = m & (df["time_min"] >= t0) & (df["time_min"] < t1)

        raw_ref = np.nanmean(df.loc[ref, "signal_mean"])
        ratio_ref = np.nanmean(df.loc[ref, "signal_bg_ratio"])

        if not np.isfinite(raw_ref):
            raw_ref = np.nanmean(df.loc[m, "signal_mean"])

        if not np.isfinite(ratio_ref) or ratio_ref == 0:
            ratio_ref = np.nanmean(df.loc[m, "signal_bg_ratio"])

        df.loc[m, "signal_corrected_raw_like"] = (
            df.loc[m, "signal_bg_ratio"] / ratio_ref * raw_ref
        )

    fig, ax = plt.subplots(figsize=(7, 5))
    fig.patch.set_facecolor("white")

    for label in df["label"].unique():
        d = df[df["label"] == label].sort_values("time_min")

        ax.plot(
            d["time_min"],
            d["signal_corrected_raw_like"],
            linewidth=2,
            label=label,
        )

    ax.axvline(0, linestyle=":", color="black", linewidth=1.2)

    ax.set_xlabel("Time since start of AMP addition (min)", fontsize=13)
    ax.set_ylabel(f"Corrected {signal_name} intensity", fontsize=13)

    title = f"Background-corrected {signal_name} signal at trench bottom"
    if AMP_title is not None:
        title = AMP_title + "\n" + title

    ax.set_title(title, fontsize=14)
    ax.legend(frameon=False, fontsize=10)
    ax.tick_params(axis="both", labelsize=12)

    plt.tight_layout()

    if save_path is not None:
        os.makedirs(os.path.dirname(save_path), exist_ok=True)
        fig.savefig(save_path, format="pdf", dpi=300, bbox_inches="tight")

    plt.show()

    return df


def plot_raw_bg_ratio_difference(
    bottom_signal_ratio_df,
    AMP_title=None,
    signal_name="LL-37",
    save_path=None,
):
    df = bottom_signal_ratio_df.copy()

    df["signal_bg_difference"] = df["signal_mean"] - df["background_mean"]

    fig, axes = plt.subplots(
        4, 1,
        figsize=(7, 10),
        sharex=True,
        constrained_layout=True,
    )

    fig.patch.set_facecolor("white")

    for label in df["label"].unique():
        df_tr = df[df["label"] == label].sort_values("time_min")

        axes[0].plot(df_tr["time_min"], df_tr["signal_mean"], linewidth=2, label=label)
        axes[1].plot(df_tr["time_min"], df_tr["background_mean"], linewidth=2, label=label)
        axes[2].plot(df_tr["time_min"], df_tr["signal_bg_ratio"], linewidth=2, label=label)
        axes[3].plot(df_tr["time_min"], df_tr["signal_bg_difference"], linewidth=2, label=label)

    for ax in axes:
        ax.axvline(0, linestyle=":", color="black", linewidth=1.2)
        ax.tick_params(axis="both", labelsize=11)

    axes[0].set_ylabel("Raw signal", fontsize=12)
    axes[1].set_ylabel("Background", fontsize=12)
    axes[2].set_ylabel("Signal / background", fontsize=12)
    axes[3].set_ylabel("Signal - background", fontsize=12)
    axes[3].set_xlabel("Time since start of AMP addition (min)", fontsize=13)

    axes[0].set_title(f"Bottom ROI {signal_name} signal", fontsize=13)
    axes[1].set_title("Side background ROI", fontsize=13)
    axes[2].set_title("Multiplicative correction", fontsize=13)
    axes[3].set_title("Additive correction", fontsize=13)

    axes[0].legend(frameon=False, fontsize=9)

    if AMP_title is not None:
        fig.suptitle(AMP_title, fontsize=14)

    if save_path is not None:
        os.makedirs(os.path.dirname(save_path), exist_ok=True)
        fig.savefig(save_path, format="pdf", dpi=300, bbox_inches="tight")

    plt.show()

    return df
