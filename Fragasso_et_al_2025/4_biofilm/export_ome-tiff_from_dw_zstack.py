# -*- coding: utf-8 -*-
"""
Bundle deconvolved per-timepoint z-stack files into corrected OME-TIFF timelapses.

For each microscopy position the script:
  1. Discovers paired GFP (c1) and SYTOX (c2) Deconwolf output files
  2. Rescales intensities by the Deconwolf calibration factor from the .log.txt file
  3. Drops the first Z_SKIP z-planes (below the glass-slide surface)
  4. Center-crops in XY to CROP_FRACTION of the original linear size
  5. Writes one TCZYX OME-TIFF per position to <in_folder>/ome_timelapse_corrected/

Output axes order: TCZYX
Channel order: [GFP, SYTOX]

@author: alessio fragasso
"""

import re
import gc
from pathlib import Path
import numpy as np
import tifffile as tiff


'''Settings'''
folder_list = [
    ## Path(r"path_to_deconwolf_output_folder_1"),
    ## Path(r"path_to_deconwolf_output_folder_2"),
]

OVERWRITE = True   # set False to skip already-written OME-TIFFs


'''Parameters'''
PX_SIZE_XY    = 0.162286        # µm/px (lateral)
PX_SIZE_Z     = 0.6             # µm/slice (axial)
CHANNEL_NAMES = ["GFP", "Sytox_Orange"]

# Z-crop: drop this many z-planes from the bottom (below glass slide)
Z_SKIP = 4

# XY center-crop: keep central CROP_FRACTION of the linear image size
# 0.5  → 50% of width and height (25% of area, 2× zoom on centre)
# 0.71 → ~50% of area
CROP_FRACTION = 0.5


'''Helpers'''
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


def build_ome_metadata():
    return {
        "axes": "TCZYX",
        "PhysicalSizeX": PX_SIZE_XY,
        "PhysicalSizeXUnit": "µm",
        "PhysicalSizeY": PX_SIZE_XY,
        "PhysicalSizeYUnit": "µm",
        "PhysicalSizeZ": PX_SIZE_Z,
        "PhysicalSizeZUnit": "µm",
        "Channel": {"Name": CHANNEL_NAMES},
    }


'''Process each folder'''
for in_folder in folder_list:

    out_folder = in_folder / "ome_timelapse_corrected"
    out_folder.mkdir(exist_ok=True)

    # Discover and group files by position and timepoint
    pattern = re.compile(r"t(?P<t>\d+)xy(?P<xy>\d+)c(?P<c>\d+)_dw\.tif")
    bundles = {}

    for f in in_folder.glob("*_dw.tif"):
        m = pattern.search(f.name)
        if m:
            c_idx = int(m.group("c"))
            if c_idx > 2:
                continue
            xy_idx = int(m.group("xy"))
            t_idx  = int(m.group("t"))
            bundles.setdefault(xy_idx, {}).setdefault(t_idx, {})[c_idx] = f

    print(f"Found {len(bundles)} positions in {in_folder.name}")

    for i, xy_idx in enumerate(sorted(bundles.keys()), start=1):

        bundle_name = f"xy{xy_idx:02d}"
        out_ome = out_folder / f"{bundle_name}_corrected.ome.tif"

        if out_ome.exists() and not OVERWRITE:
            print(f"[{i}/{len(bundles)}] Skipping (exists): {out_ome.name}")
            continue

        print(f"\n[{i}/{len(bundles)}] {bundle_name}")

        try:
            timepoints = bundles[xy_idx]
            sorted_t = sorted(timepoints.keys())
            n_t = len(sorted_t)
            n_c = 2

            # Determine dimensions from first file
            first_path = timepoints[sorted_t[0]][1]
            with tiff.TiffFile(first_path) as tf:
                z = len(tf.pages)
                y, x = tf.pages[0].shape

            if Z_SKIP >= z:
                raise ValueError(f"Z_SKIP ({Z_SKIP}) >= z-planes ({z})")

            new_z = z - Z_SKIP
            new_y = int(y * CROP_FRACTION)
            new_x = int(x * CROP_FRACTION)
            y0 = (y - new_y) // 2
            x0 = (x - new_x) // 2

            out_shape = (n_t, n_c, new_z, new_y, new_x)
            print(f"  (Z,Y,X): {z,y,x} -> cropped: {new_z,new_y,new_x} | shape: {out_shape}")

            mm = tiff.memmap(
                str(out_ome),
                shape=out_shape,
                dtype=np.uint16,
                bigtiff=True,
                ome=True,
                photometric='minisblack',
                metadata=build_ome_metadata(),
            )

            for t_step, t_idx in enumerate(sorted_t):
                for c_idx in [1, 2]:
                    f_path = timepoints[t_idx].get(c_idx)
                    if f_path is None:
                        continue

                    img = tiff.imread(str(f_path)).astype(np.float32)
                    img = img[Z_SKIP:, y0:y0 + new_y, x0:x0 + new_x]

                    sf = get_scaling_factor(f_path)
                    mm[t_step, c_idx - 1] = np.clip(img / sf, 0, 65535).astype(np.uint16)

                print(f"  wrote t{t_idx:02d} ({t_step + 1}/{n_t})", end='\r')
                gc.collect()

            mm.flush()
            del mm
            print(f"\n  Saved: {out_ome.name}")

        except Exception as e:
            print(f"\n  ERROR in {bundle_name}: {e}")

    print(f"\nDone: {in_folder.name}")
