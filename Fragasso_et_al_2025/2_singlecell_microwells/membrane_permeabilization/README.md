# 2a — Membrane-permeabilization timing and growth inhibition (CJW7845)

Estimate the timing of **outer-membrane (OM)** and **inner-membrane (IM)**
permeabilization to proteins relative to growth-rate inhibition, using strain
CJW7845 (periplasmic mScarlet-I + cytoplasmic mTagBFP2).

**Figures:** 1; 2; S2; S5; S6

## Scripts

| Script | Purpose |
|---|---|
| `get_permeabilization_time.py` | Detect OM and IM permeabilization events from the mScarlet-I and mTagBFP2 mean-intensity traces, combining each signal with its normalized first derivative. Computes the periplasm−cytoplasm mScarlet-I difference Δ(peri, cyto) to catch IM permeabilization while the OM is still intact. |
| `multi_scale_derivative.py` | Define a weighted multi-scale instantaneous growth rate (weighted average of first derivatives of Savitzky-Golay-smoothed cell area, windows 3/5/7/9/11 with weights 5/8/10/8/5). Used to mark the onset of growth inhibition and growth arrest. |

## Logic

- **OM-first phenotype:** mScarlet-I (periplasm) is lost first, then mTagBFP2
  (cytoplasm) when the IM is also permeabilized.
- **IM-first phenotype:** mScarlet-I leaks into the cytoplasm and mixes with
  mTagBFP2 → Δ(peri, cyto)_mScarlet-I goes from positive to rapidly negative.
- Growth-inhibition onset: first frame where the multi-scale growth rate drops below 0.01 min⁻¹
  and stays `< 0.012 min⁻¹` for ≥30 consecutive frames. Growth arrest: drops below 0 min⁻¹ and stays
  `< 0.005 min⁻¹` for ≥30 frames.

## Inputs

Curated per-cell feature dataFrame from `../post_processing_curation_functions.py`
(which itself consumes `../import_processing_code.py` output for CJW7845).
