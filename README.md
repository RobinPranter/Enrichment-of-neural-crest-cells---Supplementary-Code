# Enrichment of Neural Crest Cells – Supplementary Code  

This repository contains the code used to generate the results presented in the paper: “Enrichment of neural crest cells by antibody labelling and flow cytometry for single-cell transcriptomics” by Robin Pranter, Cedric Patthey, and Nathalie Feiner

## Flow Cytometry

01_plot_FACS.Rmd
Analysis and visualization of flow cytometry (FACS) data.

## RT-qPCR

qPCRAnalysis_AutoBT.Rmd
Analysis of RT-qPCR data.

## Single-Cell Transcriptomics
### 10X Chromium data

pilot3_Merged.Rmd
Reading, quality control, and filtering of 10X Chromium single-cell RNA-seq data.

### Smart-seq3 data

01_CompileData.Rmd
Compilation and initial processing of Smart-seq3 data.

02_QCFilter.Rmd
Quality control and filtering of Smart-seq3 data.

### Integrated analyses

pilot3_&SmrtSeq.Rmd
Integration of 10X Chromium and Smart-seq3 datasets, including SCT normalization, Harmony integration, and differential expression analysis.

pilot3_Merged_relExp.Rmd
Comparison of relative expression of the neural crest cell marker Sox10 across unsorted, leniently sorted, and strictly sorted wall lizard cells, as well as sorted and unsorted mouse cells from previous studies.

### Cross-species projection (scmap)

03_1st6samps_28ss_QC_scmapCao_FullData_dardel_cheap.R
Projection of strictly sorted wall lizard cells (Smart-seq3 data) onto the mouse embryonic cell atlas from Cao et al. (2019).

pilot3_scmapOnCao_FullData_dardel_cheap.R
Projection of unsorted and leniently sorted wall lizard cells (10X Chromium data) onto the mouse embryonic cell atlas from Cao et al. (2019).

scmapSolOnCao_FullData_dardel.R
Projection of sorted mouse neural crest cells from Soldatov et al. (2019) onto the mouse embryonic cell atlas from Cao et al. (2019).
