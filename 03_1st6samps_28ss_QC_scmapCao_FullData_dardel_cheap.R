
#03_1st6samps_QC_scmapCao_FullData_dardel
# Project data onto mouse atlas
#To predict NCCs in the data I will project it onto the mouse atlas from 
#Cao et al 2019 (DOI:10.1038/s41586-019-0969-x) using the program scmap 
#(https://scmap.sanger.ac.uk/scmap/)


#Set up______________________________________________________________________________

rm(list = ls()); gc()
start_time <- Sys.time()

# Load libraries
library(Seurat)
library(SingleCellExperiment)
library(magrittr) #for %>% etc
library(tidyr)    #for pivot_longer()
library(tibble)   #for rownames_to_column()
library(dplyr)    #for mutate()
library(ggplot2)
library(scmap)


#Set some global parameters:
#Number of variable features to select in scmap::selectFeatures()
#1500 seems good after testing 200-4000
features <- 1500
#Likelihood threshold for cell type assignment. 
#Should be >1/number of cell types in reference.
threshold <- 0.04


#library(future) #Already loaded with library(Seurat)
plan()
#plan("multisession", workers = 124)
plan()
options(future.globals.maxSize = 2 * 1024^4) # 2 TB available on large dardel nodes


# Define the path of the directory containing the files
path <- "/cfs/klemming/projects/snic/snic2021-23-511/WorkDirs/PmurNCC1st6_SS/"

#Source functions
source(paste0(path, "02_Analysis/pilot2_03_Functions.R"))
#____________________________________________________________________________________
print("Step 1")

## Read gene orthologies
#Read file___________________________________________________________________________
OrthoTab <-read.table(paste0(path, 
                            "01_Data/ForScmap/BM_GCRm39vsPodMur1.0_ManDL240213.txt"), 
                      sep = ",", header = TRUE, na.strings = "") %>% na.omit()
#This file was downloaded 2024-02-13 from Ensembl using BioMart with the following 
#settings:
# -Dataset:     "Mouse genes (GRCm39)"
# -Filters:     "Orthologous Common wall lizard Genes: Only"
# -Attributes:  "Common wall lizard gene name"
#               "Gene name"
#               "Gene stable ID"
#               "Common wall lizard gene stable ID"

#Rename columns
colnames(OrthoTab) <- c("PmurGeneName", "MmusGeneName", "MmusGeneStabID",
                        "PmurGeneStabID")

#Remove duplicates 
OrthoTab <- OrthoTab[!duplicated(OrthoTab$PmurGeneName), ]
OrthoTab <- OrthoTab[!duplicated(OrthoTab$MmusGeneName), ]
gc()
#____________________________________________________________________________________
print("Step 2")

## Read mouse data
### Read count table_________________________________________________________________
Mmus <- readRDS(paste0(path, "01_Data/ForScmap/gene_count_cleaned.RDS")) %>%
  as.matrix() %>% as.data.frame() #Transform sparse matrix to data frame
gc()

#Remove version number from Ensembl gene names
rownames(Mmus) <- gsub(pattern = "\\.\\d+$", 
                       replacement = "",
                       x = rownames(Mmus))
#____________________________________________________________________________________
print("Step 3")

#### Rename genes____________________________________________________________________
#Rename mouse gene stable IDs with mouse gene names
#Merge with BioMart orthology data by mouse gene stable id
Mmus$GeneStabID <- rownames(Mmus)
Mmus <- merge(x = Mmus,
              y = OrthoTab,
              by.x = "GeneStabID",
              by.y = "MmusGeneStabID")

#Set mouse gene names as rownames
rownames(Mmus) <- Mmus$MmusGeneName

#Remove uncessesary columns
Mmus <- dplyr::select(Mmus, -c("GeneStabID", "PmurGeneStabID",
                               "MmusGeneName"))
gc()
#____________________________________________________________________________________
print("Step 4")

### Read annotation file_____________________________________________________________
#Read data
Mmus_ann <- read.csv(paste0(path,
             "01_Data/ForScmap/cell_annotate.csv"))[ , c("sample", "Main_cell_type")]
gc()

#Name rows by sample
rownames(Mmus_ann) <- Mmus_ann$sample

#Subset only the cells that are in both the annotation and the count table
#Do this in both files.
Mmus_ann <- subset(Mmus_ann, sample %in% colnames(Mmus))
Mmus <- Mmus[, colnames(Mmus) %in% Mmus_ann$sample]
#ncol(Mmus) == nrow(Mmus_ann)
#ncol(Mmus)

#Remove column "sample"
Mmus_ann <- Mmus_ann["Main_cell_type"]

#Rename cell type column to be similar to scmap vignette
colnames(Mmus_ann) <- "cell_type1"
gc()
#____________________________________________________________________________________
print("Step 5")

### Create mouse sce_________________________________________________________________
#Creating the single cell experiment object for the mouse data because that is the 
#format that scmap takes.
#Inspect the reference data
#head(Mmus_ann)
#head(Mmus)[1:3, 1:3]

#Create SingleCellExperiment object
Mmus_sce <- SingleCellExperiment(assays = list(normcounts = as.matrix(Mmus)), 
                                 colData = Mmus_ann)
rm(Mmus, Mmus_ann); gc()
logcounts(Mmus_sce) <- log2(normcounts(Mmus_sce) + 1)
# use gene names as feature symbols
rowData(Mmus_sce)$feature_symbol <- rownames(Mmus_sce)
# remove features with duplicated names
Mmus_sce <- Mmus_sce[!duplicated(rownames(Mmus_sce)), ]
#Mmus_sce
#____________________________________________________________________________________
print("Step 6")

# Read and prep lizard data__________________________________________________________
#Read the data
filename <- list.files(paste0(path, "01_Data/ForScmap/"), 
                       pattern = "01_1st6samps_QC_general_output")
if (length(filename) == 1) {
  #Count table
  filt.CT <- readRDS(file = paste0(path, "01_Data/ForScmap/", filename)) %>% 
    JoinLayers(assay = "RNA") %>%
    subset(SomSt == "28ss") %>%
    LayerData(assay = "RNA", layer = "counts") %>% 
    as.data.frame()
  #Annotations
  Pmur_ann <- readRDS(file = paste0(path, "01_Data/ForScmap/", filename)
                      )[["sampleID"]] %>%as.data.frame()
} else {
  print("Choose file to read")
  print(filename)
}
gc()

print("Step 7")

#Rename lizard genes with mouse gene names
#Merge with BioMart orthology data by Podarcis gene name
filt.CT$GeneName <- rownames(filt.CT)
filt.CT <- merge(x = filt.CT,
                 y = OrthoTab,
                 by.x = "GeneName",
                 by.y = "PmurGeneName")
rm(OrthoTab); gc()
#Set mouse gene names as rownames
rownames(filt.CT) <- filt.CT$MmusGeneName

#Remove uncessesary columns
filt.CT <- dplyr::select(filt.CT, -c("GeneName", "PmurGeneStabID",
                                     "MmusGeneName", "MmusGeneStabID"))

#Subset only the cells that are in both the annotation and the count table
#Do this in both files.
Pmur_ann <- subset(Pmur_ann, rownames(Pmur_ann) %in% colnames(filt.CT))
filt.CT <- filt.CT[, colnames(filt.CT) %in% rownames(Pmur_ann)]

print("Step 8")

gc()
#____________________________________________________________________________________
print("Step 9")

### Create lizard sce________________________________________________________________
#Creating the single cell experiment object for the lizard data because that is the 
#format that scmap takes.
#Inspect the query data
filt.CT[1:3, 1:3]
head(Pmur_ann)

#Create SingleCellExperiment object
Pmur_sce <- SingleCellExperiment(assays = list(normcounts = as.matrix(filt.CT)),
                                 colData = Pmur_ann)
rm(Pmur, Pmur_ann); gc()
logcounts(Pmur_sce) <- log2(normcounts(Pmur_sce) + 1)
# use gene names as feature symbols
rowData(Pmur_sce)$feature_symbol <- rownames(Pmur_sce)
# remove features with duplicated names
Pmur_sce <- Pmur_sce[!duplicated(rownames(Pmur_sce)), ]
gc()
#____________________________________________________________________________________
print("Step 10")

## scmapCluster()
#Follow scmap vignet 
#https://bioconductor.org/packages/release/bioc/vignettes/scmap/inst/doc/scmap.html  

### Construct scmap cluster index
#### Select variable features________________________________________________________
#Select the most informative features (genes)
Mmus_sce <- selectFeatures(object = Mmus_sce, 
                           n_features = features,
                           suppress_plot = TRUE) %>% 
#table(rowData(Mmus_sce)$scmap_features)
# Construct index
#Calculate mean gene expression for every gene in each cluster
  indexCluster()
#head(metadata(Mmus_sce)$scmap_cluster_index)
#heatmap(as.matrix(metadata(Mmus_sce)$scmap_cluster_index))
#____________________________________________________________________________________
print("Step 11")

### Project lizard on mouse__________________________________________________________
#Calculate a good threshold
1/length(unique(Mmus_sce$cell_type1)) #-> threshold should be >0.027

#Project lizard cells onto mouse cluster index
Pmur_scClus_res <- scmapCluster(projection = Pmur_sce,
                                index_list = list(Cao=metadata(Mmus_sce)$
                                                    scmap_cluster_index),
                                threshold = 0.03)
#head(Pmur_scClus_res$scmap_cluster_labs)
#head(Pmur_scClus_res$scmap_cluster_siml)
#____________________________________________________________________________________
print("Step 12")

#### Inspect the results_____________________________________________________________
p1 <- SankeyForScmap(QuerryClusters = colData(Pmur_sce)$sampleID,
                     ReferenceCelltypes = Pmur_scClus_res$scmap_cluster_labs[,'Cao'],
                     Title = "1st6_SmartSeq on Cao", 
                     Subtitle = "scmapCluster()", 
                     Caption = paste0("#Lizard cells=", dim(Pmur_sce)[2], ". ",
                                      "#Mouse cells=", dim(Mmus_sce)[2], ". ", 
                                      "n_feature=1500. ",
                                      "threshold=0.03."))
ggsave(plot = p1, filename = paste0(path,"03_Results/", 
                                    "scmapCluster_1st6OnCao_28ss_FullData_dardel_",
                                    Sys.Date(), ".svg"))
#____________________________________________________________________________________
print("Step 13")

#### Rename classifications and plot again___________________________________________
Pmur_scClus_res_renamed <- Pmur_scClus_res
tempvec <- Pmur_scClus_res_renamed$scmap_cluster_labs
tempvec[tempvec=="Intermediate Mesoderm"]         <- "Not NCC"
tempvec[tempvec=="Lens"]                          <- "Not NCC"
tempvec[tempvec=="Early mesenchyme"]              <- "Possible NCC"
tempvec[tempvec=="Limb mesenchyme"]               <- "Not NCC"
tempvec[tempvec=="Neural Tube"]                   <- "Possible NCC"
tempvec[tempvec=="Schwann cell precursor"]        <- "NCC"
tempvec[tempvec=="Chondrocytes & osteoblasts"]    <- "Possible NCC"
tempvec[tempvec=="Chondroctye progenitors"]       <- "Possible NCC"
tempvec[tempvec=="Premature oligodendrocyte"]     <- "Not NCC"
tempvec[tempvec=="Cholinergic neurons"]           <- "Possible NCC"
tempvec[tempvec=="Oligodendrocyte Progenitors"]   <- "Not NCC"
tempvec[tempvec=="Cardiac muscle lineages"]       <- "Possible NCC"
tempvec[tempvec=="Jaw and tooth progenitors"]     <- "NCC"
tempvec[tempvec=="Postmitotic premature neurons"] <- "Possible NCC"
tempvec[tempvec=="Radial glia"]                   <- "Not NCC"
tempvec[tempvec=="Isthmic organizer cells"]       <- "Not NCC"
tempvec[tempvec=="Definitive erythroid lineage"]  <- "Not NCC"
tempvec[tempvec=="Sensory neurons"]               <- "Possible NCC"
tempvec[tempvec=="Notochord cells"]               <- "Not NCC"
tempvec[tempvec=="Megakaryocytes"]                <- "Not NCC"
tempvec[tempvec=="Endothelial cells"]             <- "Not NCC"
tempvec[tempvec=="White blood cells"]             <- "Not NCC"
tempvec[tempvec=="Melanocytes"]                   <- "NCC"
tempvec[tempvec=="Epithelial cells"]              <- "Not NCC"
tempvec[tempvec=="Inhibitory neurons"]            <- "Possible NCC"
tempvec[tempvec=="Neural progenitor cells"]       <- "Possible NCC"
tempvec[tempvec=="Connective tissue progenitors"] <- "Possible NCC"
tempvec[tempvec=="Ependymal cell"]                <- "Not NCC"
tempvec[tempvec=="Myocytes"]                      <- "Possible NCC"
tempvec[tempvec=="Excitatory neurons"]            <- "Possible NCC"
tempvec[tempvec=="Primitive erythroid lineage"]   <- "Not NCC"
tempvec[tempvec=="Inhibitory interneurons"]       <- "Not NCC"
tempvec[tempvec=="Stromal cells"]                 <- "Possible NCC"
tempvec[tempvec=="Inhibitory neuron progenitors"] <- "Possible NCC"
tempvec[tempvec=="Granule neurons"]               <- "Not NCC"
tempvec[tempvec=="Osteoblasts"]                   <- "Possible NCC"
tempvec[tempvec=="unassigned"]                    <- "unassigned"
tempvec[tempvec=="Hepatocytes"]                   <- "Not NCC"
Pmur_scClus_res_renamed$scmap_cluster_labs <- tempvec
rm(tempvec)

p2 <- SankeyForScmap(QuerryClusters = colData(Pmur_sce)$sampleID,
            ReferenceCelltypes = Pmur_scClus_res_renamed$scmap_cluster_labs[,'Cao'],
            Title = "1st6_SmartSeq on Cao", 
            Subtitle = "scmapCluster() renamed after proj", 
            Caption = paste0("#Lizard cells=", dim(Pmur_sce)[2], ". ",
                             "#Mouse cells=", dim(Mmus_sce)[2], ". ", 
                             "n_feature=1500. ",
                             "threshold=0.03."))
ggsave(plot = p2, filename = paste0(path,"03_Results/", 
                                    "scmapCluster_1st6OnCao_28ss_FullData_dardel_simp_",
                                    Sys.Date(), ".svg"))
#____________________________________________________________________________________
print("Step 14")

#### Count cells_____________________________________________________________________ 
#How many were classified as NCC?
length(Pmur_scClus_res_renamed$scmap_cluster_labs[Pmur_scClus_res_renamed$scmap_cluster_labs == 
                                            "NCC"])
length(Pmur_scClus_res_renamed$scmap_cluster_labs[Pmur_scClus_res_renamed$scmap_cluster_labs == 
                                            "Not NCC"])
length(Pmur_scClus_res_renamed$scmap_cluster_labs[Pmur_scClus_res_renamed$scmap_cluster_labs == 
                                            "Possible NCC"])
length(Pmur_scClus_res_renamed$scmap_cluster_labs[Pmur_scClus_res_renamed$scmap_cluster_labs == 
                                            "unassigned"])
#____________________________________________________________________________________
print("Step 15")

#Save result data____________________________________________________________________
saveRDS(list(Pmur_scClus_res, Pmur_scClus_res_renamed), paste0(path, 
                        "03_Results/",
                        "03_1st6samps_28ss_QC_scmapCao_FullData_dardel_cheap_output_", 
                        Sys.Date(),
                        ".Rds"))
#____________________________________________________________________________________

#Calculated runtime__________________________________________________________________
end_time <- Sys.time()
paste0("Elapsed time: ", end_time - start_time)
#____________________________________________________________________________________











