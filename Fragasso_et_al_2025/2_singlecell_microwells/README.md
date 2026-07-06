# 2 — Single-cell analysis in PDMS microwells

Single-cell analysis of cells immobilized in PDMS microwells and imaged by
time-lapse phase-contrast + fluorescence microscopy. This pipeline imports the
output of **OmniSegger** (drift correction + segmentation + cell linking using
the retrained Omnipose model `omnipose_phase_microwells_model`, see folder 5),
extracts per-cell features, curates trajectories, and runs three downstream
analyses: membrane-permeabilization timing, single-cell kymographs, and
intracellular DNA/ribosome perturbation.

> **The OmniSegger output is already provided in the deposited archives** — each
> `xyNN/` position folder contains the channel folders (`phase/`, `fluor1/`, …),
> `masks/`, `seg/`, the per-cell `cell/cell*.mat` files, and `clist.mat`. Start
> directly at `import_processing_code.py`; re-running OmniSegger is only needed if
> you want to re-segment from the raw images.

**Figures:** 1F, 2, 3, 4; S2–S7.

## Datasets (from S-BIAD1823)

| Dataset | Strain | Reporters | Analysis |
|---|---|---|---|
| `CJW7845_1XMIC`, `CJW7845_2XMIC`, `CJW7845_10XMIC`, `CJW7845_antibiotics` | CJW7845 | periplasmic mScarlet-I + cytoplasmic mTagBFP2 | kymographs + membrane permeabilization timing |
| `CJW7753_2XMIC` | CJW7753 | HupA-mCherry + RplA-GFP + cytoplasmic mTagBFP2 | intracellular perturbation (SCF, NC ratio) |
| `CJW7020_5-TAMRA-AMPs_AMPs_2XMIC` | CJW7020 | RplA-msfGFP | fluorescent-AMP localization (SCF*) |
| `CJW7859_5-TAMRA-AMPs_AMPs_2XMIC` | CJW7859 | HupA-msfGFP | fluorescent-AMP localization (SCF*) |
| `MG1655_5-TAMRA-AMPs_50nM_AMPs_2XMIC` | MG1655 | label-free | 5-TAMRA-AMP uptake vs. time |

## Files in this folder (shared / import / curation)

| Script | Purpose |
|---|---|
| `analysis_functions_library.py` | Shared helper functions used across the single-cell scripts. |
| `background_subtraction.py` | Local background subtraction of fluorescence images, masking dilated cell regions; GPU (cupy) accelerated. Adapted from Papagiannakis et al. |
| `import_processing_code.py` | Import linked OmniSegger per-cell MATLAB files → apply background subtraction → extract per-frame static features (area, nucleoid area, mean fluorescence, SCF, …). Outputs cropped-cell image dicts + a per-cell feature DataFrame. **Use for linked time-lapse data** (CJW7845, CJW7753, MG1655). |
| `import_processing_no_linking_code.py` | Same idea but for **un-linked** data (segmentation + drift correction only): mask labelling, per-cell cropping, and feature extraction without trajectories. **Use for CJW7020 / CJW7859** fluorescent-AMP localization (imaged every 10 min). |
| `post_processing_curation_functions.py` | Functions to compute instantaneous growth rate and curate trajectories (remove tracks that are too short, have anomalous growth rates, or lack the healthy-growth→arrest pattern). |

### Subfolders

- (membrane_permeabilization/) — timing of outer/inner
  membrane permeabilization and growth-inhibition onset (CJW7845).
- (kymographs/) — medial-axis extraction and single-cell kymographs.
- (intracellular_perturbation/) — nucleoid masking,
  nucleoid/cell length, SCF, NC ratio, and time alignment (CJW7753).

## Run order

```
[OmniSegger output provided in the archives — optional to redo:
 drift correction + segmentation + linking]
        │  per-cell MATLAB files (cell/cell*.mat, provided)
        ▼
import_processing_code.py            (linked data)
   └─ uses background_subtraction.py
        │  cropped-cell dicts + static cell feature DataFrame
        ▼
post_processing_curation_functions.py  (growth rate + trajectory curation)
        │  curated DataFrame
        ▼
  ┌─────────────────────────┬───────────────────────────┬─────────────────────────┐
  ▼                         ▼                           ▼
membrane_permeabilization/  kymographs/                 intracellular_perturbation/
(CJW7845)                   (CJW7845: Fig 2D-E, S5E,S5F) (CJW7753: Fig 3, S4)
```

For the un-linked fluorescent-AMP localization data (CJW7020 / CJW7859), use
`import_processing_no_linking_code.py` instead of the import + curation steps;
curation there is done by hard thresholds on area, focus, and fluorescence
intensity (see Methods, "Analysis of fluorescently labeled AMP uptake and
spatial localization").
