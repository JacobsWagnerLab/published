# -*- coding: utf-8 -*-
"""
Plot fluorescence kymographs from wide-trench data.

Works for any fluorescence channel (SYTOX, LL-37, etc.).
Loads outputs saved by compute_kymo_fluor_wide-trench.py and produces:
  1. Ensemble pooled overlay kymograph (signal + phase)
  2. Pooled overlay kymograph from single-trench tracks
  3. Single-trench crop + kymograph overlay panels

Requires: quantify_fluor_kymo_functions.py

@author: alessio fragasso
"""

import os
import numpy as np
import matplotlib.pyplot as plt

runfile(r"path_to_quantify_fluor_kymo_functions.py",
        wdir=r"path_to_functions_folder")


'''Settings'''
output_folder_kymo = #r"path_to_output_folder_kymo"
output_folder_fig  = #r"path_to_figure_output_folder"

master_folder_path = #r"path_to_strain_master_folder"

AMP_title  = #'AMP condition label for figure titles'
AMP_name   = #'AMP_condition_folder_name'
signal_name = #"SYTOX"   # or "LL-37", used in axis/colorbar labels

channel_name  = "sytox"   # name used when saving kymo outputs
fluor_channel = "fluor1"   # raw image subfolder: "fluor1" or "fluor2"

fr_inj = #0   # frame index at which AMP was added

os.makedirs(output_folder_fig, exist_ok=True)


'''Load ensemble kymograph outputs'''
signal_datasets = load_channel_kymo_npz_files(
    output_folder_kymo,
    channel_name=channel_name,
)

phase_datasets = load_channel_kymo_npz_files(
    output_folder_kymo,
    channel_name="phase",
)


'''Pool across positions'''
signal_frame_ids, signal_kymo, signal_2d_frames, signal_n_pos = pool_channel_datasets_by_frame(
    signal_datasets,
    ignore_zeros=True,
)

phase_frame_ids, phase_kymo, phase_2d_frames, phase_n_pos = pool_channel_datasets_by_frame(
    phase_datasets,
    ignore_zeros=True,
)

# quick diagnostic: number of positions contributing per pixel
_, _, _, _, phase_kymo_counts, _ = pool_channel_datasets_by_frame(
    phase_datasets,
    ignore_zeros=True,
    return_counts=True,
)

plt.figure(figsize=(8, 4))
plt.imshow(phase_kymo_counts.T, aspect="auto", origin="upper", cmap="viridis")
plt.colorbar(label="Number of positions contributing")
plt.title("Valid phase contributors per pixel")
plt.xlabel("Frame index")
plt.ylabel("y position")
plt.show()

print("Frame IDs match:", np.array_equal(signal_frame_ids, phase_frame_ids))

trench_rel_box = signal_datasets[0]["metadata"]["trench_rel_box"]
print("trench_rel_box:", trench_rel_box)


'''Plot ensemble overlay kymograph'''
pdf_path = os.path.join(
    output_folder_fig,
    f"{AMP_name}_{signal_name}_phase_overlay_pooled_kymograph.pdf",
)

plot_pooled_overlay_kymograph(
    pooled_frame_ids=signal_frame_ids,
    pooled_kymo_signal=signal_kymo,
    pooled_kymo_phase=phase_kymo,
    fr_inj=fr_inj,
    trench_rel_box=trench_rel_box,
    AMP_title=AMP_title,
    signal_name=signal_name,

    overlay_mode="colormap",   # "green" or "colormap"
    overlay_cmap="inferno",

    manual_signal_vmin=None,   # set to override percentile-based limits
    manual_signal_vmax=None,
    signal_vmin_prc=1,
    signal_vmax_prc=90,

    manual_phase_vmin=None,
    manual_phase_vmax=None,
    phase_vmin_prc=0,
    phase_vmax_prc=95,

    phase_invert=False,
    signal_alpha=0.85,
    zero_background="light_grey",
    show_trench_lines=False,

    ylims_px=(250, 1200),
    ytick_step_um=10,
    save_path=pdf_path,
)


'''Load single-trench tracks'''
single_tracks_signal = load_single_trench_tracks(
    output_folder_kymo,
    channel_name=channel_name,
)

single_tracks_phase = load_single_trench_tracks(
    output_folder_kymo,
    channel_name="phase",
)


'''Pool single-trench tracks and plot overlay kymograph'''
signal_frame_ids_st, signal_kymo_st, signal_counts_st, chosen_signal = pool_single_tracks_by_frame(
    single_tracks_signal,
    ignore_zeros=True,
    min_valid_tracks=3,
)

phase_frame_ids_st, phase_kymo_st, phase_counts_st, chosen_phase = pool_single_tracks_by_frame(
    single_tracks_phase,
    ignore_zeros=True,
    min_valid_tracks=3,
)

print("Frame IDs match (single tracks):", np.array_equal(signal_frame_ids_st, phase_frame_ids_st))
print("Signal kymo shape:", signal_kymo_st.shape, "| Phase kymo shape:", phase_kymo_st.shape)

pdf_path_singletrack_pool = os.path.join(
    output_folder_fig,
    f"{AMP_name}_{signal_name}_phase_overlay_pooled_from_single_tracks.pdf",
)

plot_pooled_overlay_kymograph(
    pooled_frame_ids=signal_frame_ids_st,
    pooled_kymo_signal=signal_kymo_st,
    pooled_kymo_phase=phase_kymo_st,
    fr_inj=fr_inj,
    trench_rel_box=single_tracks_signal[0]["metadata"]["trench_rel_box"],
    AMP_title=AMP_title + "\npooled from single-trench tracks",
    signal_name=signal_name,

    overlay_mode="colormap",
    overlay_cmap="inferno",

    manual_signal_vmin=None,
    manual_signal_vmax=None,
    signal_vmin_prc=1,
    signal_vmax_prc=90,

    manual_phase_vmin=None,
    manual_phase_vmax=None,
    phase_vmin_prc=0,
    phase_vmax_prc=95,

    phase_invert=False,
    signal_alpha=0.85,
    zero_background="light_grey",
    show_trench_lines=False,

    ylims_px=(250, 1300),
    ytick_step_um=10,
    save_path=pdf_path_singletrack_pool,
)


'''Extract tracked trench crops'''
frames_to_extract = [
    fr_inj,         # t = 0
    fr_inj + 30,    # t = +30
    fr_inj + 60,    # t = +60
    fr_inj + 100,   # t = +100
    ## fr_inj + 150,
]

track_idx_list = list(range(len(single_tracks_signal)))

crop_records_dict = {}

for track_index in track_idx_list:
    crop_records_dict[track_index] = extract_tracked_trench_crops_for_frames(
        single_tracks_sytox=single_tracks_signal,
        track_index=track_index,
        frames_to_extract=frames_to_extract,
        output_folder_kymo=output_folder_kymo,
        master_folder_path=master_folder_path,
        fluor_channel=fluor_channel,
        pad_x=75,
        pad_top=250,
        pad_bottom=100,
        max_frame_delta=2,
        use_sorted=False,
        signal_bkg_sub=False,   # match saved kymo setting
    )

    print("Processed crops for track", track_index)


'''Plot single-trench crop + kymograph overlays'''
phase_lookup = make_phase_track_lookup(single_tracks_phase)

pos_to_plot_list = [
    ## "xy01",
    ## "xy02",
    ## "xy03",
]

for pos_to_plot in pos_to_plot_list:

    track_indices = find_track_indices_by_pos(single_tracks_signal, pos_to_plot)
    print(f"Track indices for {pos_to_plot}:", track_indices)

    for track_idx in track_indices:

        tr_s = single_tracks_signal[track_idx]
        key = get_track_match_key(tr_s)

        if key not in phase_lookup:
            print("No matching phase track for:", key)
            continue

        phase_idx = phase_lookup[key]
        crop_records_this_track = crop_records_dict[track_idx]

        if len(crop_records_this_track) == 0:
            print("No crops for track", track_idx)
            continue

        track_id = int(tr_s["track_id"])
        pos = tr_s["metadata"].get("pos")

        for crop_i in range(len(crop_records_this_track)):

            frame_used = int(crop_records_this_track[crop_i]["frame"])

            pdf_path_single = os.path.join(
                output_folder_fig,
                f"{AMP_name}_{pos}_track_{track_id:03d}_frame_{frame_used:04d}_crop_kymo_overlay.pdf",
            )

            plot_selected_crop_and_track_kymo(
                crop_records=crop_records_this_track,
                single_tracks_sytox=[tr_s],
                single_tracks_phase=[single_tracks_phase[phase_idx]],
                crop_record_index=crop_i,
                track_index=0,
                fr_inj=fr_inj,
                AMP_title=AMP_title + f", track = {track_idx}, crop = {crop_i}",
                signal_name=signal_name,

                manual_signal_vmin=None,   # set to fix contrast across tracks
                manual_signal_vmax=None,
                signal_vmin_prc=1,
                signal_vmax_prc=95,

                crop_phase_vmin_prc=1,
                crop_phase_vmax_prc=99.5,
                kymo_phase_vmin_prc=1,
                kymo_phase_vmax_prc=99.5,
                manual_crop_phase_vmin=None,
                manual_crop_phase_vmax=None,
                manual_kymo_phase_vmin=None,
                manual_kymo_phase_vmax=None,

                phase_invert=False,
                signal_alpha=0.9,
                zero_background="light_grey",

                overlay_mode="colormap",
                overlay_cmap="inferno",
                scale_bar_um=5,

                kymo_tlims_min=(0, 200),
                ylims_px=(50, 1000),
                ytick_step_um=10,

                save_path=pdf_path_single,
                show=True,
            )


