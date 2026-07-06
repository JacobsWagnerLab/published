# 3 — Wide-trench microfluidics

Analysis of bacteria imaged in wide-trench (mother-machine-style) microfluidic
devices. Two independent readouts are produced from the same drift-corrected
image sequences:

- **Growth rate** — *single-cell* analysis from segmented, linked cells
  (growth-rate kymographs).
- **SYTOX Green / LL-37 fluorescence** — *not single-cell*; fluorescence
  kymographs quantified at the trench level.

Segmentation/linking (growth-rate workflow only) is done with **OmniSegger**
using the retrained wide-trench Omnipose models (`omnipose_ribo_wide-trench_model`
for growth rate; `omnipose_phase_wide-trench_model` to define cell-containing
regions when building phase / (SYTOX/LL-37 fluor) overlay kymographs— see folder 5).

> **The OmniSegger output is already provided in the deposited archives** (per-
> position channel folders, `masks/`, `seg/`, per-cell `cell/cell*.mat`, and
> `clist.mat`). Start at `import_processing_code_wide-trench.py`; re-running
> OmniSegger is only needed if you want to re-segment from the raw images.

**Figures:** growth rate (CJW7753) → 5C–F, 5I, S8A–B; SYTOX Green (CJW7845) →
5G, S8C–D; LL-37 (CJW7859) → 5H, S8E–G.
**Videos:** 5–6 (growth rate); 7-8 (overlay SYTOX Green/phase); 9 (overlay LL-37/phase)

## Datasets (from S-BIAD1823)

| Dataset | Strain | Channels used | Analysis |
|---|---|---|---|
| `CJW7753_2XMIC_wide-trench`, `CJW7753_10XMIC_wide-trench` | CJW7753 | ribosome (RplA-GFP) + phase | growth-rate kymographs ([`growth_rate/`](growth_rate/)) |
| `CJW7845_2XMIC_wide-trench` | CJW7845 | phase + SYTOX Green | SYTOX Green fluorescence kymographs ([`fluorescence_kymographs/`](fluorescence_kymographs/)) |
| `CJW7859_4xMIC_wide-trench_LL-37_fluor` | CJW7859 | phase + LL-37 | LL-37 fluorescence kymographs + bottom-of-trench quantification ([`fluorescence_kymographs/`](fluorescence_kymographs/)) |

## Files in this folder

| Script | Purpose |
|---|---|
| `drift_correction_wide_trenches.py` | Correct slow residual y-drift in the wide-trench phase-contrast image sequences (per-trench ROI selection) after SuperSegger drift correction was run. Shared first step for both workflows. |
| `import_processing_code_wide-trench.py` | (Growth rate) Import the linked OmniSegger per-cell MATLAB files and build per-cell DataFrames. |
| `post_processing_curation_code_wide-trench.py` | (Growth rate) Curate trajectories: remove tracks < 15 min, with anomalous normalized growth rate (> 0.04 or < −0.03 min⁻¹), or large center-of-mass jumps (> 50 px). |


### Subfolders

- [`growth_rate/`](growth_rate/) — growth-rate kymographs, single-trench
  plots/videos, and mp4 export.
- [`fluorescence_kymographs/`](fluorescence_kymographs/) — SYTOX Green / LL-37
  fluorescence kymographs and the LL-37 bottom-of-trench quantification.

## Run order

Both workflows start from the drift-corrected image sequences and then diverge:

```
drift_correction_wide_trenches.py        (shared: residual y-drift correction)
        │
        ├─────────────────────────────────────────────┐
        ▼                                               ▼
  GROWTH RATE (single-cell)                    FLUORESCENCE (not single-cell)
        │                                               │
  [OmniSegger output provided —                 compute_kymo_fluor_wide-trench.py
   optional to redo:                            + remaining fluorescence_kymographs/
   segmentation + linking]                        functions
        ▼
  import_processing_code_wide-trench.py
        ▼
  post_processing_curation_code_wide-trench.py
        ▼
  growth_rate/
```

Time is offset relative to AMP injection: **t = 0 = AMP addition** throughout.
