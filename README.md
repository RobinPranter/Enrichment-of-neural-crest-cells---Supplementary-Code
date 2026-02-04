# Enrichment of Neural Crest Cells – Supplementary Code  

This repository contains the code and metadata and some data used to generate the results presented in the paper: “Enrichment of neural crest cells by antibody labelling and flow cytometry for single-cell transcriptomics” by Robin Pranter, Cedric Patthey, and Nathalie Feiner. 
The scripts assume that you are standing in a directory named "~/WorkDirs/PmurNCCPilot_10X_rerun/" which contains the three directories 01_Data/, 02_Analysis/ and 03_Results.


## Flow Cytometry

FlowCytometry.Rmd
Analysis and visualization of flow cytometry (FACS) data. Reads the data and metadata in 01_Data/FlowCyt and produces the plots in figure 3. 

## RT-qPCR

qPCRAnalysis_AutoBT.Rmd
Analysis of RT-qPCR data. Reads the data in 01_Data/ and produces the plots in figure 4 and results in table 1.


## Single-Cell Transcriptomics
01_Data/ should contain the metadata files provided here, and the count matrices found at [insert GEO accession]. For the final three scripts, you further need to acquire data from previous studies (Soldatov 2019 DIO:10.1126/science.aas9536; and Cao 2019 DIO:10.1038/s41586-019-0969-x). 
02_Analysis/ should contain the scripts provided here. 
03_Results is where figures etc will be saved. 
 
The scripts are numbered in the order they should be run. Bellow is a flow chart of the pipeline followed by a short description of each script. 

```text
01_CompileSmrtS.Rmd    03_CompileAndFilter10X.Rmd
          |               |
          v               |
02_FilterSmrtS.Rmd        |
          |               |
          +-------+-------+
                  |
                  v
04_IntAnalysis10XAndSmrtS.Rmd
05_RelExp10XAndSmrtS.Rmd
```

01_CompileSmrtS.Rmd - Compilation and initial processing of Smart-seq3 data. 

02_FilterSmrtS.Rmd - Quality control and filtering of Smart-seq3 data. 

03_CompileAndFilter10X.Rmd - Compilation, quality control, and filtering of 10X Chromium single-cell RNA-seq data. 

04_IntAnalysis10XAndSmrtS.Rmd -  Integrated analysis of 10X Chromium and Smart-seq3 datasets, including SCT normalization, Harmony integration, and differential expression analysis. 

05_RelExp10XAndSmrtS.Rmd - Comparison of relative expression of the neural crest cell marker Sox10 across unsorted, leniently sorted, and strictly sorted wall lizard cells, as well as sorted and unsorted mouse cells from previous studies. 

### Cross species projection (scmap)  
Three scripts that each use the scmap method from Kiselev et al 2018 (DOI:10.1038/nmeth.4644) to projects scRNA-seq data onto the whole mouse embryo cell atlas from Cao et al. (2019).  
  
06_scmap10XOnCao.R - Projection of unsorted and leniently sorted wall lizard cells (10X Chromium data) onto the mouse embryonic cell atlas from Cao et al. (2019).  
  
07_scmapSolOnCao.R - Projection of sorted mouse neural crest cells from Soldatov et al. (2019) onto the mouse embryonic cell atlas from Cao et al. (2019).  
  
08_scmapSmrtSOnCao.R -  Projection of strictly sorted wall lizard cells (Smart-seq3 data) onto the mouse embryonic cell atlas from Cao et al. (2019).  

