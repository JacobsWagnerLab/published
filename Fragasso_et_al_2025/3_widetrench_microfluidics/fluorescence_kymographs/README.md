# 3b — Wide-trench SYTOX Green / LL-37 fluorescence kymographs

Fluorescence kymographs of SYTOX Green and LL-37 along
the wide trenches, and the quantification of when LL-37 fluorescence reaches the
bottom of empty trenches.

**Figures:** 5G–H; S8C–G.

## Scripts

| Script | Purpose |
|---|---|
| `quantify_fluor_kymo_functions.py` | **Shared functions library** for the fluorescence-kymograph pipeline: trench detection, ROI extraction, single-trench tracking, ensemble/single-trench kymograph computation, save/load, pooling, overlay plotting, and LL-37 bottom-of-trench signal helpers. Imported by the scripts below. |
| `compute_kymo_fluor_wide-trench.py` | Background-subtract fluorescence + phase; restrict phase to cell-containing regions using dilated Omnipose masks; detect trench edges/top/bottom from the phase x-projection; project signals onto the trench y-axis (excluding 50 px per side); build and save single-trench (tracked) and ensemble kymographs. |
| `plot_kymo_fluor.py` | Load saved kymographs, pool across positions, and plot fluorescence kymographs overlaid on matched phase kymographs (ensemble and single-trench). Generic over channel (SYTOX Green or LL-37). |
| `quantify_LL-37_fluor_empty_trench.py` | Quantify mean LL-37 fluorescence in a 40-px band at the bottom of selected tracked **empty** trenches, background-corrected with a matched region outside the trench, aligned to AMP injection (Fig. S8F–G). |

## Notes

- Single-trench kymographs link trenches across frames by x-center position and
  keep tracks present for ≥20 frames; ensemble kymographs average y-profiles
  across trenches/positions per frame.
- Time is aligned to AMP injection (**t = 0 = AMP addition**).
- The LL-37 dataset was saved under the channel name `"sytox"` in the original
  scripts; `plot_kymo_fluor.py` therefore loads it via that same name (set
  `channel_name`/`signal_name` at the top of the script accordingly).

## Run order

```
compute_kymo_fluor_wide-trench.py   →  saved kymographs (per position)
        ▼
plot_kymo_fluor.py                  →  overlay kymograph figures
quantify_LL-37_fluor_empty_trench.py →  bottom-of-trench LL-37 traces
            (both import quantify_fluor_kymo_functions.py)
```
