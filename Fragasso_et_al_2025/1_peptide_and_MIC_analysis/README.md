# 1 — Peptide physicochemical properties & MIC / growth-rate analysis

Code for the bulk-level characterization of AMPs/PMs: their
physicochemical properties and their minimum inhibitory concentrations (MICs)
measured by microplate reader.

**Figures:** 1A–B; S1A–C, F–G; S9A–B; S10B.

## Scripts

| Script | Purpose |
|---|---|
| `protein_analysis.py` | Main analysis script for the peptide panel. Defines the peptides (sequences, source, structure, amidation, disulfide bridges, MIC, mechanism of action), builds the per-peptide property table via `get_protein_df`, and produces both figures and tables: the AMP summary table (**Fig. 1A**), the property scatterplots — MW vs. net charge, Lys vs. Arg fraction, GRAVY vs. proline fraction, hydrophobic vs. polar fraction (**Fig. 1B**), the physicochemical-properties table (**Fig. S1A**), and the sequences / class-and-type table (**Fig. S10B**). |
| `protein_analysis_functions.py` | Function library imported by `protein_analysis.py`. Defines `get_protein_df(label, aa_seq, amidated, disulfide_bridges, source)`, which computes GRAVY index (Kyte-Doolittle), polar/hydrophobic/cationic fractions, net charge (K, R, C-terminal amidation), and molecular weight (peptide-bond, amidation, and disulfide-bridge corrections). |
| `microplate_reader_analysis.py` | From microplate-reader OD₆₀₀ curves, compute the instantaneous growth rate (linear fit of log(OD) over a 12-point / 60-min moving window), the maximum growth rate, and the MIC. |

## Data

- **Peptide analysis:** amino-acid sequences are defined inside
  `protein_analysis.py` (no external data file needed).
- **MIC / growth rate:** microplate-reader OD₆₀₀ time-series exported as text
  files (sampled every 5 min, up to 48 h), for MG1655 (or ATCC 25922 for the
  biofilm strain MIC), mixed 50:50 with titrated AMP concentrations. The MIC is
  the lowest AMP concentration giving a maximum growth rate `< 0.2 log(OD)/h`
  within the first 20 h; the **mode** across ≥3 replicates is the final MIC.

## Scripts can be run independently

1. **Peptide properties:** run `protein_analysis.py`. It imports
   `protein_analysis_functions.py` (same folder). By default figures are only
   displayed; set `save_path` at the top of the script to also write PNG + PDF.
2. **MIC / growth rate:** run `microplate_reader_analysis.py` on the exported OD
   text files (see `plot_max_gr_rate_vs_conc` and `get_max_gr_params`).
