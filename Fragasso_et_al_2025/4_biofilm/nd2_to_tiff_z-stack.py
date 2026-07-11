# -*- coding: utf-8 -*-
"""
Export per-(T,P,C) TIFFs from ND2 without loading large stacks into RAM.
Memory-optimized version: resets internal SDK buffers and recreates Dask views 
per-loop to prevent "hidden" memory accumulation.
"""

import argparse
import gc
import sys
from pathlib import Path

import numpy as np
import tifffile as tiff
import nd2


def main():
    ap = argparse.ArgumentParser(
        description="Export TIFFs from ND2 splitting T/P/C while streaming planes."
    )

    ap.add_argument("nd2_path", type=Path, help="Path to input .nd2")
    ap.add_argument("-o", "--out", type=Path, required=True, help="Output directory")

    zgroup = ap.add_mutually_exclusive_group()
    zgroup.add_argument("--z-index", type=int, default=None)
    zgroup.add_argument("--z-max", action="store_true")
    zgroup.add_argument("--z-stack", action="store_true")

    ap.add_argument("--start-t", type=int, default=0)
    ap.add_argument("--end-t", type=int, default=None)
    ap.add_argument("--start-p", type=int, default=0)
    ap.add_argument("--end-p", type=int, default=None)
    ap.add_argument("--start-c", type=int, default=0)
    ap.add_argument("--end-c", type=int, default=None)
    ap.add_argument("--c-one-based", action="store_true")
    ap.add_argument("--overwrite", action="store_true")
    ap.add_argument("--bigtiff", action="store_true")

    args = ap.parse_args()
    args.out.mkdir(parents=True, exist_ok=True)
    base = args.nd2_path.stem

    with nd2.ND2File(args.nd2_path) as f:
        # metadata extraction
        axes = list(f.sizes.keys())
        sizes = dict(f.sizes)
        dtype = f.dtype

        print(f"ND2 axes: {axes} sizes: {sizes}", flush=True)

        T, P, C, Z = sizes.get("T", 1), sizes.get("P", 1), sizes.get("C", 1), sizes.get("Z", 1)
        Y, X = sizes.get("Y"), sizes.get("X")

        t_end = args.end_t if "T" in sizes and args.end_t is not None else T
        p_end = args.end_p if "P" in sizes and args.end_p is not None else P
        c_end = args.end_c if "C" in sizes and args.end_c is not None else C

        # Determine Z mode
        if "Z" in sizes:
            if args.z_stack: z_mode = "stack"
            elif args.z_max: z_mode = "max"
            else:
                z_mode = "index"
                z_index_val = args.z_index if args.z_index is not None else 0
        else:
            z_mode, z_index_val = "none", None

        # --- Helper Functions ---
        def get_plane(ti, pi, ci, zi, arr_view):
            idx = []
            for ax in axes:
                if ax == "T": idx.append(ti)
                elif ax == "P": idx.append(pi)
                elif ax == "C": idx.append(ci)
                elif ax == "Z": idx.append(zi if zi is not None else 0)
                elif ax in ("Y", "X"): idx.append(slice(None))
                else: idx.append(0)
            
            # Compute and force a hard copy to detach from the memory map
            plane = np.array(arr_view[tuple(idx)].compute(scheduler="single-threaded")).copy()
            
            # Grayscale conversion if RGB
            if plane.ndim == 3 and plane.shape[-1] in (3, 4):
                plane = (0.2126*plane[...,0] + 0.7152*plane[...,1] + 0.0722*plane[...,2]).astype(dtype)
            
            return np.squeeze(plane)

        def write_stack(path, ti, pi, ci, arr_view):
            tmp = path.with_suffix(".part.tif")
            is_big = args.bigtiff or (int(Z)*int(Y)*int(X)*np.dtype(dtype).itemsize > 3.5e9)
            with tiff.TiffWriter(tmp, bigtiff=is_big) as tw:
                for zi in range(Z):
                    tw.write(get_plane(ti, pi, ci, zi, arr_view), photometric="minisblack", contiguous=True)
            tmp.replace(path)

        def write_max(path, ti, pi, ci, arr_view):
            proj = None
            for zi in range(Z):
                plane = get_plane(ti, pi, ci, zi, arr_view)
                proj = np.maximum(proj, plane) if proj is not None else plane.copy()
            tiff.imwrite(path, proj, photometric="minisblack")

        # --- Main Loop ---
        for ti in range(args.start_t, t_end):
            for pi in range(args.start_p, p_end):
                for ci in range(args.start_c, c_end):
                    
                    # 1. Create a FRESH Dask view for this specific TPC 
                    # This prevents the Dask graph from accumulating across iterations
                    current_arr = f.to_dask()

                    t_lab = f"t{ti+1:02d}"
                    p_lab = f"xy{pi+1:02d}"
                    c_val = ci+1 if args.c_one_based else ci
                    out_path = args.out / f"{base}_{t_lab}{p_lab}c{c_val}.tif"

                    if out_path.exists() and not args.overwrite:
                        continue

                    print(f"Exporting {out_path.name}...", flush=True)

                    if z_mode == "stack":
                        write_stack(out_path, ti, pi, ci, current_arr)
                    elif z_mode == "max":
                        write_max(out_path, ti, pi, ci, current_arr)
                    else:
                        zi = z_index_val if z_mode == "index" else None
                        tiff.imwrite(out_path, get_plane(ti, pi, ci, zi, current_arr), photometric="minisblack")

                    # 2. MANDATORY CLEANUP
                    del current_arr
                    gc.collect()
                    
                    # 3. NUKING THE SDK CACHE
                    # This reaches into the nd2 library to clear C-level buffers
                    if hasattr(f, '_rdr'):
                        try:
                            # Re-initializing the reader handle forces buffer release
                            f._rdr.close()
                            f._rdr.open()
                        except:
                            pass 

    print("\nAll exports finished. Memory cleared.")

if __name__ == "__main__":
    main()