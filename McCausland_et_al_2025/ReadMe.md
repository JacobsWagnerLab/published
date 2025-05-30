# Code for McCausland et al. 2025

The code here contains relevant analysis files and code necessary for generating main figures in [McCausland et al., 2025 on bioRxiv](https://www.biorxiv.org/content/10.1101/2025.01.08.631998v1?ct=). Here ther are three folders and 1 analysis file. 

- Borrelia_Cell_Segmentation, authored by Joshua W. McCausland. This is an analysis package to segment Borrelia cells from phase contrast images, including tools for scanning the medial axis of segemented objects for demograph generations and basic shape description. This relies primarily on Otsu segmentation to identify cells. This analysis generated the HADA cell dataframe used for **Figures 3D-F**.
- Muropeptide_library, authored by Irnov Irnov. This generates the theoretical muropeptide library used throughout all the mass spec analysis in the paper. 
- XCMS_screen_supernatants, written primarily by Irnov Irnov with some modification from Joshua W. McCausland. This is a basic R script to analyze mzML files from supernatant data and identify peaks of interest.
- figure_generation, written by Irnov Irnov (**Figure 1A**) and Joshua W. McCausland (everything except **Figure 1A**). Small datasets are included in this folder (no mass spec or imaging experiments). All mass spec and imaging data will be uploaded to GlycoPOST and Biostudies, respectively.
  - Main figure MS data: https://glycopost.glycosmos.org/entry/GPST000537
  - Revision experiment MS data (extra patients in Figure 6 and S14, dead cell PG experiment in Figure S13): https://glycopost.glycosmos.org/entry/GPST000595
  - Raw microscopy data for HADA imaging in Figure 3: https://www.ebi.ac.uk/biostudies/bioimages/studies/S-BIAD1573
