# 4 — Biofilm 3D analysis

Quantification and 3D rendering of biofilm widefield fluorescence z-stacks for
strain **ATCC 25922-GFP** (constitutive GFP), with SYTOX Orange (dead-cell
marker) or EbbaBiolight 680 (amyloid/matrix optotracer).

Z-stacks are first deconvolved with **Deconwolf** (Wernersson et al., 2024;
https://deconwolf.fht.org/), 50 iterations, GPU; theoretical PSFs generated at
xy = 162.286 nm, z = 600 nm, 145 z-slices, NA = 0.95, n = 1, emission 525 nm
(GFP) / 610 nm (SYTOX Orange) / 640 nm (EbbaBiolight 680). The Deconwolf scaling
factors in the `.log.txt` files are reapplied before quantification.

The exact theoretical PSFs used for deconvolution are provided at S-BIAD1823 in
`PSF_theoretical.zip` — one Deconwolf-generated `.tiff` per channel
(`PSF_GFP_0.6um.tiff`, `PSF_SYTOX_Orange_0.6um.tiff`,
`PSF_EbbaBioLight680_0.6um.tiff`), each with its `.log.txt` recording the
generation parameters above.

**Figures:** 6A–C; S9C–D. **Videos:** biofilm 3D montages (Videos 10-17)

## Datasets (from S-BIAD1823)

`ATCC25922GFP_biofilm_rep1`, `ATCC25922GFP_biofilm_rep2`, `ATCC25922GFP_biofilm_rep3`
— three biological replicates of ATCC 25922-GFP biofilms (control + AMP-treated
wells), each containing the raw widefield z-stacks (per timepoint, per
channel). Deconvolve with Deconwolf before running the scripts below.

The per-position → condition mapping is listed in each archive's
`notes_<date>...txt` file (2 positions per well). AMPs are added ~1 h into
acquisition (after 2 frames), together with SYTOX Orange (50 nM). Positions
**xy01–xy14** are the control + AMP wells imaged with GFP + SYTOX Orange;
**xy15–xy16** are the control wells imaged with GFP + EbbaBiolight 680 (matrix
optotracer) — the same split used by `export_3D_imaris_montage.py`.

`PSF_theoretical.zip` (also at S-BIAD1823) contains the per-channel theoretical
PSFs used to deconvolve these z-stacks (see the deconvolution parameters above).

## Scripts

| Script | Purpose |
|---|---|
| `quantify_biofilm_stats.py` | From paired deconvolved GFP (c1) + SYTOX (c2) `*_dw.tif` stacks: rescale by Deconwolf factor, drop the first 3 z-planes and keep the next 40, segment GFP biomass by adaptive local threshold and SYTOX by fixed threshold (although not used in the final quantification plots), and compute per-timepoint biovolume, mean intensity, integrated and total raw signal. Saves per-position pickle DataFrames. |
| `export_ome-tiff_from_dw_zstack.py` | Bundle per-timepoint deconvolved z-stacks into 5D **OME-TIFF** timelapses (TCZYX; channels GFP + SYTOX Orange). Applies Deconwolf scaling, stores physical voxel size, drops the first 4 z-planes, and keeps the central 50% of the field of view in x/y. |
| `export_3D_imaris_montage.py` | Assemble Imaris-exported TIFF render sequences into montage mp4 videos (side + top views, scale bars, colorbars, timestamps). **xy01–xy14 = GFP + SYTOX Orange** (3-panel); **xy15–xy16 = GFP + EbbaBiolight 680** (4-panel). |

## Run order

```
Widefield z-stacks (per timepoint, per channel)
        │  Deconwolf (3D deconvolution, external; theoretical PSFs in PSF_theoretical.zip)
        ▼  *_dw.tif + *.log.txt
   ┌──────────────────────────────┬──────────────────────────────────┐
   ▼                              ▼
quantify_biofilm_stats.py        export_ome-tiff_from_dw_zstack.py
(biomass / Fig 6C, S9D)           │  5D OME-TIFF
                                  ▼
                          Imaris (.ims convert → 3D render → TIFF sequences, external)
                                  ▼
                          export_3D_imaris_montage.py  → montage videos
```

## Notes

- `export_3D_imaris_montage.py` expects one subfolder per rendered channel/view
  inside its `ROOT_DIR`, named with suffixes `_gfp`, `_sytox`, `_comp_top`
  (SYTOX positions) or `_gfp`, `_ebba_glow`, `_gfp_top`, `_ebba_glow_top`
  (EBBA positions). For the EbbaBiolight 680 channel it uses a Fiji-style Glow
  LUT; place a `glow.lut` next to the script to use the exact LUT, otherwise a
  built-in approximation is used.
- Edit the `path_to_...` / `ROOT_DIR` placeholders at the top of each script
  before running.
