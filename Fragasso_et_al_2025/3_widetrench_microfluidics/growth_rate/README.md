# 3a — Wide-trench growth-rate kymographs

Single-cell growth-rate kymographs along the wide trenches,
plus single-trench plots and mp4 videos.

**Figures:** 5C–F, 5I; S8A–B. **Videos:** 5–6.

## Scripts

| Script | Purpose |
|---|---|
| `compute_kymo_growth_rate_wide-trench.py` | From curated single-cell trajectories, estimate instantaneous growth rates, bin the cell center-of-mass y-coordinate into 12-px (0.79 µm) bins, prune isolated bins with an iterative neighborhood filter, and build mean-normalized-growth-rate kymographs over time × position. Pools across trenches/replicates. |
| `plot_single_trench_growth_rate.py` | Combine single-cell growth-rate data from an individual trench with the corresponding image sequence to produce single-trench figures/videos (Fig. 5C). |
| `export_mp4_wide-trench.py` | Encode the rendered frame sequences into mp4 videos (Videos 5–6), via `imageio_ffmpeg`. |

## Inputs

Curated wide-trench DataFrame from
`../post_processing_curation_code_wide-trench.py` and the cropped trench image
sequences from `../drift_correction_wide_trenches.py`.
