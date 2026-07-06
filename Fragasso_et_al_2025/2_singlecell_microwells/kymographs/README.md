# 2b — Single-cell kymographs

Compute single-cell kymographs by projecting fluorescence signals onto the
cell medial axis.

**Figures:** 2D, 2E, S5E, S5F; S3.

## Scripts

| Order | Script | Purpose |
|---|---|---|
| 1 | `get_medial_axis.py` | Compute the medial axis of each cell mask: 20× upsample, smooth, `skeletonize`, find poles by linear regression of the skeleton extremities (first/last 3/8), then fit a 3rd-order polynomial through skeleton + poles. |
| 2 | `get_1D_proj_medial.py` | Map each cell pixel (x, y) to medial-axis coordinates (l = arc length, d = distance from axis). Returns the mapped dataFrames and the mean-intensity projection along the axis (keeping pixels with d < 4 px). |
| 3 | `kymograph_single_cell.py` | Bin the sub-pixel medial-axis length into original-pixel bins (0.065 µm/px), average pixels with the same l, and stack the resulting 1D projections over time to form the kymograph. |

## Inputs

Cropped-cell masks and fluorescence channels from the import step
(`../import_processing_code.py`), for CJW7753 (HupA-mCherry / RplA-GFP) and
related strains.
