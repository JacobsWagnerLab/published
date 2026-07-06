# -*- coding: utf-8 -*-
"""
Quantify GFP biomass and SYTOX Green signal from deconvolved wide-field biofilm z-stacks.

For each microscopy position the script:
  1. Reads paired Deconwolf-deconvolved GFP (c1) and SYTOX (c2) z-stacks
  2. Rescales intensities by the Deconwolf calibration factor from the .log.txt file
  3. Truncates the z-stack to a fixed region above the glass slide
  4. Segments GFP biomass with adaptive local thresholding; SYTOX with a fixed threshold
  5. Computes per-timepoint biovolume, mean intensity, integrated signal, and raw volume signal
  6. Saves per-position results as pickle DataFrames under <dw_folder>/biomass_extracted_dfs/

Output DataFrame columns (one row per timepoint per position):
  file_base, position_id, timepoint, biomass_volume_gfp, biomass_volume_sytox,
  mean_gfp, mean_sytox, integrated_biomass_gfp, integrated_biomass_sytox,
  gfp_initial_object_count, total_raw_signal_gfp, total_raw_signal_sytox

@author: alessio fragasso
"""

import re
import tifffile as tiff
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
from pathlib import Path
from scipy import ndimage as ndi
from skimage.morphology import remove_small_objects
from skimage.filters import threshold_local


'''Settings'''
target_directory = #r"path_to_deconwolf_output_folder"   # folder containing *_dw.tif files


'''Parameters'''
PX_SIZE_XY = 0.162286   # µm/px (lateral)
PX_SIZE_Z  = 0.6        # µm/slice (axial)
VOXEL_VOL  = PX_SIZE_XY * PX_SIZE_XY * PX_SIZE_Z  # µm³

# z-stack truncation: skip first Z_START_INDEX slices (below glass), keep Z_SLICES_COUNT
Z_START_INDEX  = 3
Z_SLICES_COUNT = 40
Z_END_INDEX    = Z_START_INDEX + Z_SLICES_COUNT

# GFP segmentation: adaptive local threshold
GFP_SIGMA      = 10.0
GFP_BLOCK_SIZE = 61
GFP_OFFSET     = 80
GFP_FLOOR      = 300   # minimum threshold value (retains expanding biofilm edges)

# SYTOX segmentation: fixed threshold
SYTOX_SIGMA  = 1.5
SYTOX_THRESH = 500

# Minimum object size [voxels] for both channels
SAFE_MIN_SIZE_PX = 200

# Maximum timepoints to process (0-indexed; 48 = 24 h at 30-min intervals)
MAX_TIMEPOINTS = 48

# Diagnostic monitoring z-slice (raw index before truncation)
MONITOR_Z_RAW = 10


'''Helper: parse Deconwolf scaling factor from log file'''
def get_scaling_factor(tif_path):
    log_path = Path(str(tif_path) + ".log.txt")
    if not log_path.exists():
        return 1.0
    try:
        with open(log_path, 'r') as f:
            matches = re.findall(r"scaling:\s+([0-9.]+)", f.read())
            if matches:
                return float(matches[-1])
    except Exception:
        pass
    return 1.0


'''Main processing function'''
def batch_process_by_position(dw_folder_path):
    """
    Iterate over all positions found in dw_folder_path, process each timepoint,
    and save per-position result DataFrames to <dw_folder>/biomass_extracted_dfs/.
    """

    dw_folder = Path(dw_folder_path)
    output_subfolder = dw_folder / "biomass_extracted_dfs"
    output_subfolder.mkdir(exist_ok=True)

    all_c1_files = sorted(dw_folder.glob("*_t*xy*c1_dw.tif"))

    if not all_c1_files:
        print("No Deconwolf Channel 1 (*c1_dw.tif) files found.")
        return

    # Group GFP files by position
    position_bundles = {}
    for f_path in all_c1_files:
        m = re.search(r'_t(\d+)xy(\d+)c1_dw\.tif', f_path.name)
        if m:
            pos_id = f"xy{m.group(2)}"
            position_bundles.setdefault(pos_id, []).append(f_path)

    print(f"Found {len(position_bundles)} positions.")
    print(f"Output: {output_subfolder}\n")

    for pos_id, c1_files in sorted(position_bundles.items()):

        output_pkl = output_subfolder / f"biomass_results_{pos_id}.pkl"
        if output_pkl.exists():
            print(f"Skipping {pos_id} (already processed)")
            continue

        print(f"=== Processing {pos_id} ===")
        position_results = []

        for gfp_path in sorted(c1_files):

            m = re.search(r'_t(\d+)xy\d+c1_dw\.tif', gfp_path.name)
            if not m:
                continue

            raw_t   = int(m.group(1))
            timepoint = raw_t - 1   # 0-indexed: _t01 → timepoint 0

            if timepoint > MAX_TIMEPOINTS:
                print(f"  Reached time limit ({MAX_TIMEPOINTS} frames). Stopping.")
                break

            sytox_path = gfp_path.with_name(
                gfp_path.name.replace("c1_dw.tif", "c2_dw.tif")
            )
            if not sytox_path.exists():
                print(f"  Missing SYTOX file for t={raw_t}. Skipping.")
                continue

            # Read and rescale
            gfp_raw   = tiff.imread(gfp_path).astype(np.float32)
            sytox_raw = tiff.imread(sytox_path).astype(np.float32)

            sf_gfp   = get_scaling_factor(gfp_path)
            sf_sytox = get_scaling_factor(sytox_path)

            gfp_3d   = np.clip(gfp_raw   / sf_gfp,   0, 65535).astype(np.uint16)[Z_START_INDEX:Z_END_INDEX]
            sytox_3d = np.clip(sytox_raw / sf_sytox,  0, 65535).astype(np.uint16)[Z_START_INDEX:Z_END_INDEX]

            # Unmasked volume signal (backup metric)
            total_raw_gfp   = float(np.sum(gfp_3d))
            total_raw_sytox = float(np.sum(sytox_3d))

            # Smooth
            gfp_filt   = ndi.gaussian_filter(gfp_3d,   sigma=GFP_SIGMA)
            sytox_filt = ndi.gaussian_filter(sytox_3d, sigma=SYTOX_SIGMA)

            # Segment GFP: adaptive local threshold
            local_thresh = threshold_local(
                gfp_filt,
                block_size=GFP_BLOCK_SIZE,
                method='gaussian',
                offset=GFP_OFFSET,
            )
            local_thresh = np.clip(local_thresh, a_min=GFP_FLOOR, a_max=800)
            mask_gfp = remove_small_objects(
                gfp_filt >= local_thresh,
                min_size=SAFE_MIN_SIZE_PX,
            )

            # Segment SYTOX: fixed threshold
            mask_sytox = remove_small_objects(
                sytox_filt >= SYTOX_THRESH,
                min_size=SAFE_MIN_SIZE_PX,
            )

            # Label GFP objects for object count
            _, n_objects_gfp = ndi.label(mask_gfp)

            # Biovolume and mean intensity
            n_vox_gfp   = int(np.sum(mask_gfp))
            n_vox_sytox = int(np.sum(mask_sytox))

            vol_gfp   = n_vox_gfp   * VOXEL_VOL
            vol_sytox = n_vox_sytox * VOXEL_VOL

            mean_gfp   = float(np.mean(gfp_3d[mask_gfp]))     if n_vox_gfp   > 0 else 0.0
            mean_sytox = float(np.mean(sytox_3d[mask_sytox])) if n_vox_sytox > 0 else 0.0

            # Diagnostic plot at monitor z-slice
            plot_z = MONITOR_Z_RAW - Z_START_INDEX
            if plot_z < 0 or plot_z >= gfp_3d.shape[0]:
                plot_z = gfp_3d.shape[0] // 2

            fig, axes = plt.subplots(2, 3, figsize=(19, 10))
            fig.suptitle(
                f"{pos_id} | t={timepoint:02d} ({timepoint * 0.5:.1f} h) | z={plot_z + Z_START_INDEX}",
                fontsize=13, fontweight='bold',
            )

            axes[0, 0].imshow(gfp_3d[plot_z],           cmap='gray',   vmin=10,  vmax=1000)
            axes[0, 0].set_title("GFP raw (corrected)")
            axes[0, 1].imshow(mask_gfp[plot_z],          cmap='viridis')
            axes[0, 1].set_title(f"GFP mask | {n_objects_gfp} objects")
            axes[0, 2].imshow(mask_gfp[plot_z] * gfp_3d[plot_z], cmap='gray', vmin=10, vmax=1000)
            axes[0, 2].set_title("GFP masked overlay")

            axes[1, 0].imshow(sytox_3d[plot_z],          cmap='gray',   vmin=100, vmax=3000)
            axes[1, 0].set_title("SYTOX raw (corrected)")
            axes[1, 1].imshow(mask_sytox[plot_z],         cmap='plasma')
            axes[1, 1].set_title(f"SYTOX mask | thresh={SYTOX_THRESH}")
            axes[1, 2].imshow(mask_sytox[plot_z] * sytox_3d[plot_z], cmap='gray', vmin=100, vmax=3000)
            axes[1, 2].set_title("SYTOX masked overlay")

            for ax in axes.ravel():
                ax.axis('off')

            plt.tight_layout()
            plt.show()

            position_results.append({
                "file_base":               gfp_path.name.split("_t")[0],
                "position_id":             pos_id,
                "timepoint":               timepoint,
                "biomass_volume_gfp":      vol_gfp,
                "biomass_volume_sytox":    vol_sytox,
                "mean_gfp":                mean_gfp,
                "mean_sytox":              mean_sytox,
                "integrated_biomass_gfp":  mean_gfp   * vol_gfp,
                "integrated_biomass_sytox":mean_sytox * vol_sytox,
                "gfp_initial_object_count":n_objects_gfp,
                "total_raw_signal_gfp":    total_raw_gfp,
                "total_raw_signal_sytox":  total_raw_sytox,
            })

            print(
                f"  t={timepoint:02d} ({timepoint * 0.5:.1f} h) | "
                f"Vol GFP: {vol_gfp:.1f} µm³ | Raw GFP: {total_raw_gfp:.2e}"
            )

        if position_results:
            pd.DataFrame(position_results).to_pickle(output_pkl)
            print(f"--> Saved: {output_pkl.name}\n")
        else:
            print(f"--> No results for {pos_id}.\n")

    print("Done.")


'''Run'''
batch_process_by_position(target_directory)
