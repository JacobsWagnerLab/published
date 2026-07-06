# Fragasso et al. — Antimicrobial peptide (AMP) single-cell, wide-trench microfluidics, and biofilm imaging analysis

Analysis code accompanying the manuscript titled 'Time-resolved phenotyping at 
subcellular resolution reveals shared principles and key trade-offs across 
antimicrobial peptide activities'

This repository contains the custom Python code used to process and quantify
the microscopy data. The raw microscopy images are deposited
separately on the **BioImage Archive** under accession
**[S-BIAD1823](https://www.ebi.ac.uk/biostudies/bioimages/studies/S-BIAD1823)**.

---

## 1. What is here

The code is grouped into five pipelines, each in its own numbered folder with a
dedicated `README.md`:

| Folder | Pipeline | Main figures |
|--------|----------|--------------|
| [`1_peptide_and_MIC_analysis/`](1_peptide_and_MIC_analysis/) | Peptide physicochemical properties (Fig. S1A) + sequences/classification (Fig. S10B) + MIC / max growth-rate from microplate-reader OD curves | Fig. 1A–B; S1A–C, F–G; S9A–B; S10B |
| [`2_singlecell_microwells/`](2_singlecell_microwells/) | Single-cell analysis in PDMS microwells: import of OmniSegger output, background subtraction, curation, membrane-permeabilization timing, single-cell kymographs, intracellular DNA/ribosome perturbation | Fig. 1F, 2–4; S2–S7 |
| [`3_widetrench_microfluidics/`](3_widetrench_microfluidics/) | Wide-trench microfluidics: drift correction, cropping, growth-rate kymographs, and SYTOX Green / LL-37 fluorescence kymographs | Fig. 5; S8 |
| [`4_biofilm/`](4_biofilm/) | Biofilm widefield z-stacks: Deconwolf-deconvolved biomass quantification, OME-TIFF assembly, and Imaris 3D montage videos | Fig. 6; S9C–D |
| [`5_omnipose_retraining/`](5_omnipose_retraining/) | Generation of curated training data and the retrained Omnipose segmentation models | (segmentation for all single-cell data) |

---

## 2. Data: from BioImage Archive to analysis

All imaging data live at BioImage Archive accession
**[S-BIAD1823](https://www.ebi.ac.uk/biostudies/bioimages/studies/S-BIAD1823)**.
Datasets are organized by *E. coli* strain and treatment. Each dataset folder
contains per-condition `.zip` archives (one archive per replicate / AMP /
concentration).

> **Each archive includes a `notes_<date>[_<strain>].txt` file** at its top
> level, recording the acquisition metadata for that condition: date, strain,
> temperature, growth medium, frame interval, AMP (or antibiotic) and
> concentration, imaged positions (`xyNN`), and the **frame at which the AMP was
> added / injected** (needed to set t = 0 in the analysis). Consult this file
> before running the scripts to fill in condition-specific parameters (e.g. the
> AMP-addition frame).

> **OmniSegger output is already included — you do not need to run OmniSegger.**
> For all single-cell and wide-trench datasets, each archive unzips to one folder
> per imaged position (`xyNN/`) that already contains the **complete OmniSegger
> output**: the raw channel folders (`phase/`, `fluor1/`, `fluor2/`, `fluor3/`),
> the segmentation `masks/`, the intermediate `seg/` files, the per-cell linked
> MATLAB files (`cell/cell*.mat`), and the position cell list (`clist.mat`).
> You can therefore start directly at the Python **import** step (folder 2 / 3).
> Re-running OmniSegger with the retrained models (folder 5) is only necessary if
> you want to re-segment from scratch.

| BioStudies dataset folder | Strain | Reporters | Used by pipeline | Figures |
|---|---|---|---|---|
| `CJW7845_1XMIC`, `CJW7845_2XMIC`, `CJW7845_10XMIC`, `CJW7845_antibiotics` | CJW7845 | periplasmic mScarlet-I + cytoplasmic mTagBFP2 | 2 — membrane permeabilization | Fig. 1F, 2; S2, S6 |
| `CJW7753_2XMIC` | CJW7753 | HupA-mCherry + RplA-GFP + cytoplasmic mTagBFP2 | 2 — kymographs + intracellular perturbation | Fig. 2D–E, 3, S5; S3, S4 |
| `CJW7020_5-TAMRA-AMPs_AMPs_2XMIC` | CJW7020 | RplA-msfGFP | 2 — fluorescent-AMP localization (no linking) | Fig. 4C; S7 |
| `CJW7859_5-TAMRA-AMPs_AMPs_2XMIC` | CJW7859 | HupA-msfGFP | 2 — fluorescent-AMP localization (no linking) | Fig. 4C; S7 |
| `MG1655_5-TAMRA-AMPs_50nM_AMPs_2XMIC` | MG1655 | label-free (5-TAMRA-AMP spike-in) | 2 — AMP uptake over time (with linking) | Fig. 4A; S7A |
| `CJW7753_2XMIC_wide-trench`, `CJW7753_10XMIC_wide-trench` | CJW7753 | ribosome (RplA-GFP) + phase | 3 — growth-rate kymographs | Fig. 5C–F, I; S8A–B |
| `CJW7845_2XMIC_wide-trench` | CJW7845 | phase + SYTOX Green | 3 — SYTOX Green fluorescence kymographs | Fig. 5G; S8C–D |
| `CJW7859_4xMIC_wide-trench_LL-37_fluor` | CJW7859 | phase + LL-37 | 3 — LL-37 fluorescence kymographs | Fig. 5H; S8E–G |
| `ATCC25922GFP_biofilm_rep1`, `ATCC25922GFP_biofilm_rep2`, `ATCC25922GFP_biofilm_rep3` | ATCC 25922-GFP | GFP + SYTOX Orange / EbbaBiolight 680 | 4 — biofilm z-stacks | Fig. 6; S9C–D |
| `PSF_theoretical` | — | theoretical PSFs (GFP / SYTOX Orange / EbbaBiolight 680) | 4 — Deconwolf deconvolution of the biofilm z-stacks | Fig. 6; S9C–D |

### General data → analysis flow

```
BioImage Archive (.zip per condition)
        │  download + unzip
        ▼
Per-position folders (xyNN/) with channel images + masks + per-cell MATLAB files
        │  (already the OmniSegger output — see note above)
        │
        ├─ Single-cell & wide-trench ───────────────────────────────┐
        │     [OmniSegger already run — optional to redo:            │
        │      drift correction, segmentation with retrained         │
        │      Omnipose models (folder 5), and cell linking]         │
        │            │  per-cell MATLAB files (provided)             │
        │            ▼                                               │
        │     import  →  background subtraction  →  feature          │
        │     extraction (folder 2 / 3)                              │
        │            │  pickled DataFrames + cropped-cell dicts      │
        │            ▼                                               │
        │     curation  →  quantification / kymographs / plots       │
        │                                                            │
        └─ Biofilm ──────────────────────────────────────────────────┘
              Deconwolf (3D deconvolution)
                     │  deconvolved z-stacks (*_dw.tif + .log.txt)
                     ▼
              export_ome-tiff_from_dw_zstack.py  /  quantify_biofilm_stats.py
                     │
                     ▼
              Imaris (3D render → TIFF sequences) → export_3D_imaris_montage.py
```

See each folder's `README.md` for the exact script order, inputs/outputs, and
parameters for that pipeline.

---

## 3. Upstream tools (not in this repository)

These tools run before / alongside the Python code and must be installed
separately. **For single-cell and wide-trench data, OmniSegger and Omnipose are
optional** — their output is already included in the archived datasets (see the
note in Section 2), so they are only needed if you want to re-segment from the
raw images. Deconwolf and Imaris are required to reproduce the biofilm pipeline
from the deposited z-stacks.

- **OmniSegger** (drift correction, segmentation, cell linking) —
  https://github.com/tlo-bot/omnisegger (MATLAB; combines SuperSegger +
  Omnipose). *Optional — output provided.*
- **Omnipose** (segmentation backend) —
  https://github.com/kevinjohncutler/omnipose . Use with the retrained models
  in [`5_omnipose_retraining/retrained_omnipose_models/`](5_omnipose_retraining/).
  *Optional — only for re-segmenting from raw images.*
- **Deconwolf** (biofilm 3D deconvolution) — https://deconwolf.fht.org/
  (Wernersson et al., 2024).
- **Imaris 10.2 + ImarisFileConverter 11.0.1** (biofilm 3D rendering, commercial).

---

## 4. Software environment

All custom code was developed in **Python 3.9**. Key packages:

```
numpy        pandas        scipy         matplotlib    seaborn
scikit-image tifffile      imageio       imageio-ffmpeg
Pillow       keyboard	   cupy*         torch*        pims (ND2Reader_SDK)* 
```

`*` GPU-accelerated / hardware- or reader-specific:
- **cupy** (+ a CUDA toolkit) accelerates large-array operations
  (`background_subtraction.py`, medial-axis computation). An NVIDIA A40 GPU was
  used in this study; on a CPU-only machine, replace `cupy`/`cupyx` calls with
  their `numpy`/`scipy` equivalents.
- **torch** is used during the import step.
- **pims** with the Nikon ND2 SDK reads `.nd2` files.

There is no single environment file; install the packages above into a fresh
conda/virtualenv (Python 3.9). Scripts are run individually (e.g. in Spyder /
an IDE), not as a single command-line pipeline — most expose a `Settings`
block of paths and parameters near the top that you edit before running.

---

## 5. How to run, in short

1. Download the relevant dataset(s) from
   [S-BIAD1823](https://www.ebi.ac.uk/biostudies/bioimages/studies/S-BIAD1823)
   and unzip. For single-cell / wide-trench data, each archive already contains
   the full OmniSegger output per position (see Section 2).
2. Open the relevant pipeline folder, read its `README.md`, edit the path
   placeholders (`path_to_...`) at the top of each script, and run the scripts
   in the documented order — starting at the **import** step, which reads the
   provided per-cell MATLAB files directly.
3. *(Optional)* To re-segment from the raw images instead of using the provided
   masks/cell files, run the data through **OmniSegger** with the matching
   retrained Omnipose model (folder 5) before the import step.
4. For biofilm data, deconvolve with **Deconwolf** first, then run the folder 4
   scripts.


---

## 6. Citation & data availability

- **Code:** the Jacobs-Wagner lab repository,
  https://github.com/JacobsWagnerLab/published/tree/master/Fragasso_et_al_2025.
- **Imaging data:** BioImage Archive accession
  [S-BIAD1823](https://www.ebi.ac.uk/biostudies/bioimages/studies/S-BIAD1823).

Please cite the manuscript and the BioImage Archive accession when using this
code or data.
