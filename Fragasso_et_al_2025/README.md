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
**[S-BIAD1823](https://www.ebi.ac.uk/biostudies/bioimages/studies/S-BIAD1823)** —
the naming convention and a strain → pipeline → figure mapping are below.

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

The data are deposited as a **flat list of per-condition `.zip` archives** (not a
nested folder tree), each named
`<strain>_<AMP>_<concentration>[_rep<N>][_wide-trench].zip`
(e.g. `CJW7753_CecA_1uM_rep1.zip`). Concentrations in the filenames are the
**absolute** values used. The per-AMP **1× MIC** values are given in **Fig. 1A**;
the 2×, 4×, and 10× MIC conditions are simply those values multiplied by the
corresponding factor. Search the S-BIAD1823 file list by strain + AMP. The table
below maps each **strain / experiment type** to its reporters, the pipeline that
uses it, and the figures.

| BioStudies archive name pattern | Strain | Reporters | Used by pipeline | Figures |
|---|---|---|---|---|
| `CJW7845_<AMP>_<conc>[_rep<N>]` | CJW7845 | periplasmic mScarlet-I + cytoplasmic mTagBFP2 | 2 — membrane permeabilization, single-cell kymographs | Fig. 1F, 2; S2, S3, S5, S6 |
| `CJW7753_<AMP>_<conc>[_rep<N>]` | CJW7753 | HupA-mCherry + RplA-GFP + cytoplasmic mTagBFP2 | 2 — intracellular perturbation | Fig. 3, S4 |
| `CJW7020_<AMP>_<conc>[_rep<N>]` (5-TAMRA-AMPs) | CJW7020 | RplA-msfGFP | 2 — fluorescent-AMP localization (no linking) | Fig. 4C; S7 |
| `CJW7859_<AMP>_<conc>[_rep<N>]` (5-TAMRA-AMPs) | CJW7859 | HupA-msfGFP | 2 — fluorescent-AMP localization (no linking) | Fig. 4C; S7 |
| `MG1655_<AMP>_<conc>[_rep<N>]` (5-TAMRA-AMPs) | MG1655 | label-free (5-TAMRA-AMP spike-in) | 2 — AMP uptake over time (with linking) | Fig. 4A; S7A |
| `CJW7753_<AMP>_<conc>_wide-trench` | CJW7753 | ribosome (RplA-GFP) + phase | 3 — growth-rate kymographs | Fig. 5C–F, I; S8A–B |
| `CJW7845_<AMP>_<conc>_wide-trench` | CJW7845 | phase + SYTOX Green | 3 — SYTOX Green fluorescence kymographs | Fig. 5G; S8C–D |
| `CJW7859_<AMP>_<conc>_wide-trench` (LL-37) | CJW7859 | phase + LL-37 (intrinsic fluorescence) | 3 — LL-37 fluorescence kymographs | Fig. 5H; S8E–G |
| `ATCC25922GFP_biofilm_rep<N>` | ATCC 25922-GFP | GFP + SYTOX Orange / EbbaBiolight 680 | 4 — biofilm quantification + 3D rendering | Fig. 6; S9C–D |
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

The code was developed and tested with **Python 3.10** (3.10.16) on **Windows 10/11**,
using an NVIDIA GPU with **CUDA 12.1** (an NVIDIA A40 was used in this study). The
GPU-accelerated steps additionally require CUDA-capable hardware (see notes below).

### Tested package versions

| Package | Version | Used for |
|---|---|---|
| python | 3.10.16 | — |
| numpy | 2.2.6 | all pipelines |
| pandas | 2.2.3 | all pipelines |
| scipy | 1.15.2 | all pipelines |
| matplotlib | 3.10.1 | plotting (needs ≥3.10 for the `berlin` colormap) |
| seaborn | 0.13.2 | plotting |
| scikit-image | 0.25.2 | segmentation / morphology |
| tifffile | 2025.3.13 | TIFF / OME-TIFF I/O (biofilm) |
| imageio | 2.37.0 | image reading |
| imageio-ffmpeg | 0.6.0 | mp4 video export |
| pillow | 10.4.0 | image I/O |
| keyboard | 0.13.5 | interactive curation |
| cupy-cuda12x | 13.4.0 | GPU array ops (needs CUDA 12.x) |
| torch | 2.5.1 (cu121) | import step |
| torchvision | 0.20.1 | with torch |
| pims + pims-nd2 | 0.7 / 1.1 | reading `.nd2` (needs Nikon ND2 SDK) |
| nd2 | 0.10.4 | SDK-free `.nd2` reading (biofilm ND2 → TIFF) |

`*` Hardware- / reader-specific notes:
- **cupy** (+ a CUDA 12.x toolkit) accelerates large-array operations
  (`background_subtraction.py`, medial-axis computation). On a CPU-only machine,
  replace `cupy`/`cupyx` calls with their `numpy`/`scipy` equivalents.
- **torch** is used during the import step.
- **pims + pims-nd2** require the proprietary **Nikon ND2 SDK** (a native library
  installed separately) to read `.nd2` files. The pure-Python **nd2** package reads
  `.nd2` without the SDK and is used for the biofilm ND2 → TIFF conversion.

### Installing the environment

```bash
# 1. Create and activate a fresh environment
conda create -n amp_env python=3.10
conda activate amp_env

# 2. Install PyTorch (CUDA 12.1 build)
pip install torch==2.5.1 torchvision==0.20.1 --index-url https://download.pytorch.org/whl/cu121

# 3. Install the remaining pinned packages
pip install -r requirements.txt
```

Typical install time on a normal desktop is ~10–15 min (dominated by the PyTorch
download). `cupy-cuda12x` requires a working CUDA 12.x installation; on a CPU-only
machine, omit it from `requirements.txt` and adapt the `cupy`/`cupyx` calls as noted
above. `pims-nd2` will import only if the Nikon ND2 SDK is installed — it is needed
only for reading raw `.nd2` files, not for the deposited (already-extracted) data.

Scripts are run individually (e.g. in Spyder / Jupyter / an IDE), not as a single
command-line pipeline — most expose a `Settings` block of paths and parameters near
the top that you edit before running. See [`6_demos/`](6_demos/) for runnable,
self-contained examples on small datasets.

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

### Runnable demos

To make each pipeline easy to try, [`6_demos/`](6_demos/) provides a runnable
demo for every pipeline, packaged as a Jupyter notebook that walks through the
steps end-to-end and shows the expected output inline. The **"Data needed"**
column indicates whether a demo runs from the provided intermediate dataframes
alone, or additionally needs the raw images:

| Demo | Pipeline it demonstrates | Data needed |
|---|---|---|
| [`1_import_curation_single-cell_analysis`](6_demos/1_import_curation_single-cell_analysis/) | Single-cell microwells: import, background subtraction, curation | **Raw** (import reads OmniSegger output) |
| [`2_import_curation_single-cell_analysis_no_linking`](6_demos/2_import_curation_single-cell_analysis_no_linking/) | Un-linked single-cell feature extraction (fluorescent-AMP localization) | **Raw** (import reads OmniSegger output) |
| [`3_import_curation_single-cell_analysis_wide-trenches`](6_demos/3_import_curation_single-cell_analysis_wide-trenches/) | Wide-trench single-cell import + curation | **Raw** (import reads OmniSegger output) |
| [`4_CJW7845_growth_and_permeabilization_kinetics_medial_axis_kymograph`](6_demos/4_CJW7845_growth_and_permeabilization_kinetics_medial_axis_kymograph/) | Growth-arrest / membrane-permeabilization timing + medial-axis kymographs | Dataframes only |
| [`5_CJW7753_growth_inhibition_and_intracellular_perturbations`](6_demos/5_CJW7753_growth_inhibition_and_intracellular_perturbations/) | Nucleoid/ribosome perturbation, NC ratio, cell length | Dataframes only |
| [`6_CJW7020_colocalization_fluor_AMPs`](6_demos/6_CJW7020_colocalization_fluor_AMPs/) | Fluorescent-AMP colocalization / SCF | Dataframes only |
| [`7_CJW7753_single-cell_growth_rate_analysis_wide-trench`](6_demos/7_CJW7753_single-cell_growth_rate_analysis_wide-trench/) | Wide-trench single-cell growth-rate kymographs | Dataframes only |
| [`8_CJW7845_quantify_sytox_fluor_kymo_wide-trench`](6_demos/8_CJW7845_quantify_sytox_fluor_kymo_wide-trench/) | Wide-trench SYTOX / fluorescence kymographs | Dataframes for the pooled/ensemble kymographs; **raw** for the single-trench crop panels |
| [`9_ATCC25922GFP_biofilm_deconvolution_and_quantification`](6_demos/9_ATCC25922GFP_biofilm_deconvolution_and_quantification/) | Biofilm ND2 → TIFF, Deconwolf deconvolution, biomass quantification | Dataframes for the biomass plots; **raw** z-stacks for the ND2→TIFF / deconvolution / quantification steps |

**Getting the demo data.** Two sources, depending on how far back you want to start:

1. **Intermediate dataframes (`demos_datasets.zip`)** — download from BioStudies
   [S-BIAD1823](https://www.ebi.ac.uk/biostudies/bioimages/studies/S-BIAD1823) and
   unzip it **at the repository root**; it merges into `6_demos/`, placing each
   demo's `output*/` dataframes next to its notebook. This is all that **demos 4–7**
   need to run top to bottom, and it drives the pooled kymographs (demo 8) and the
   biomass plots (demo 9).
2. **Raw images** — needed for the **import** step of **demos 1–3** (and for the
   raw-image steps of demos 8 and 9). Download the corresponding per-condition
   archive from [S-BIAD1823](https://www.ebi.ac.uk/biostudies/bioimages/studies/S-BIAD1823)
   (e.g. `CJW7753_CecA_1uM_rep1.zip`) and unzip it into that demo's folder so the
   **experiment folder** (e.g. `CJW7753_CecA_1uM_rep1/`, which holds the `xy*/`
   position folders and the `notes_*.txt`) sits next to the notebook — matching
   the `experiment_path` set at the top of the notebook. The notebook's early
   markdown cells note which archive to fetch.

Every notebook also ships with its outputs saved inline, so you can see the
expected results (figures and printed values) **without running anything** — the
data downloads above are only needed if you want to re-execute a demo.


---

## 6. Citation & data availability

- **Code:** the Jacobs-Wagner lab repository,
  https://github.com/JacobsWagnerLab/published/tree/master/Fragasso_et_al_2025.
- **Imaging data:** BioImage Archive accession
  [S-BIAD1823](https://www.ebi.ac.uk/biostudies/bioimages/studies/S-BIAD1823).

Please cite the manuscript and the BioImage Archive accession when using this
code or data.
