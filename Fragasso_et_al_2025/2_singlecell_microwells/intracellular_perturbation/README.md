# 2c — Intracellular DNA/ribosome perturbation (CJW7753)

Quantify AMP-induced intracellular reorganization of DNA and ribosomes using
strain CJW7753 (HupA-mCherry + RplA-GFP): the signal correlation factor (SCF),
the nucleocytoplasmic (NC) ratio, and time alignment of trajectories by the SCF
peak.

**Figures:** 3C, 3D, 3F, 3G; S4.

## Scripts

| Script | Purpose |
|---|---|
| `nucleoid_mask.py` | Compute the nucleoid mask from the HupA (DNA) fluorescence (`get_fluor_mask`, Otsu-based thresholding). |
| `nucleoid_cell_length.py` | Project nucleoid pixels onto the cell medial axis and compute the arc length (= nucleoid length); the same approach on the full medial axis gives the cell length, hence the NC ratio. |
| `time_alignment.py` | Build new time axes by aligning trajectories to peak values of chosen markers (`area`, `SCF`, NC ratio) within a time window. Supports the multi-step alignment used for class I vs. class II peptides (align by peak area or min NC ratio, then by max SCF). |

## Key quantities

- **SCF**: pixel-based Spearman correlation between RplA-GFP and HupA-mCherry
  over a 4-px eroded cell mask (computed during the import step,
  `../import_processing_code.py`).
- **NC ratio**: nucleoid length / cell length.
- **Alignment**: trajectories are synchronized so that peak SCF is at t = 0
  (Methods, "Analysis of intracellular perturbations upon AMP uptake").

## Inputs

Curated per-cell DataFrame (`../post_processing_curation_functions.py`) with SCF
already computed at import; nucleoid masks/lengths are added here, then
trajectories are time-aligned for the ensemble plots.
