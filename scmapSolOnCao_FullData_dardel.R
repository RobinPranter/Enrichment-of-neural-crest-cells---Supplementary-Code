
# Set up_____________________________________________________________________________
#start_time <- Sys.time()
rm(list = ls()); invisible(gc())
path = "/cfs/klemming/projects/snic/snic2021-23-511/WorkDirs/PmurNCCPilot_10X/"

library(magrittr)
library(SingleCellExperiment)
library(scmap)
library(ggplot2)

options(future.globals.maxSize = 2 * 1024^4) # 2 TB available on large dardel nodes

source(paste0(path, "02_Analysis/pilot2_03_Functions.R"))
#____________________________________________________________________________________


# Read gene orthology table__________________________________________________________

#It may seem counter intuitive to read the lizard-mouse orthology table since I will only 
#use mouse data. But I will use the mouse column to arrange the count tables to have 
#the same rownames.

#Read file
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
print("Step 1")

# Read Cao data 
## Count table_______________________________________________________________________
Cao <- readRDS(paste0(path, "01_Data/ForScmap/gene_count_cleaned.RDS")) %>%
  as.matrix() %>% as.data.frame() #Transform sparse matrix to data frame
gc()

#Transform sparse matrix to data frame
#Further subsample by randomly selecting a fraction of all cells to make it workable
#Using all 100000 cells the script takes ~30 min so it is not too bad. 
#Cao <- as.data.frame(as.matrix(Cao#[ , round(runif(ncol(Cao))-0.30) == 1]))
#invisible(gc())

#Remove version number from Ensembl gene names
rownames(Cao) <- gsub(pattern = "\\.\\d+$", 
                      replacement = "",
                      x = rownames(Cao))
#____________________________________________________________________________________
print("Step 2")

### Rename genes_____________________________________________________________________
#Rename mouse gene stable IDs with mouse gene names
#Merge with BioMart orthology data by mouse gene stable id
Cao$GeneStabID <- rownames(Cao)
Cao <- merge(x = Cao,
             y = OrthoTab,
             by.x = "GeneStabID",
             by.y = "MmusGeneStabID")

#Set mouse gene names as rownames
rownames(Cao) <- Cao$MmusGeneName

#Remove uncessesary columns
Cao <- dplyr::select(Cao, -c("GeneStabID", "PmurGeneStabID","MmusGeneName"))
gc()
#____________________________________________________________________________________
print("Step 3")

## Annotations_______________________________________________________________________
#Read data
Cao_ann <- read.csv(paste0(path, "01_Data/ForScmap/", 
                           "cell_annotate.csv"))[ , c("sample", "Main_cell_type")]
gc()

#Name rows by sample
rownames(Cao_ann) <- Cao_ann$sample

#Subset only the cells that are in both the annotation and the count table
#Do this in both files.
Cao_ann <- subset(Cao_ann, sample %in% colnames(Cao))
Cao <- Cao[, colnames(Cao) %in% Cao_ann$sample]
#ncol(Cao) == nrow(Cao_ann)
#ncol(Cao)

#Remove column "sample"
Cao_ann <- Cao_ann["Main_cell_type"]

#Rename cell type column to be similar to scmap vignette
colnames(Cao_ann) <- "cell_type1"
gc()
#____________________________________________________________________________________
print("step 4")

# Read Soldatov data
## Count table#______________________________________________________________________
#Read count tables
Sol_1 <- read.table(paste0(path, "01_Data/ForScmap/", 
                           "E10.5_post-otic_Sox10_counts.txt"), 
                    sep = " ")
Sol_1$gene <- rownames(Sol_1)
Sol_2 <- read.table(paste0(path, "01_Data/ForScmap/", 
                           "E10.5_tail_Wnt1_counts.txt"), 
                    sep = " ")
Sol_2$gene <- rownames(Sol_2)
Sol_3 <- read.table(paste0(path, "01_Data/ForScmap/", 
                           "E10.5_trunk_Wnt1_counts.txt"), 
                    sep = " ")
Sol_3$gene <- rownames(Sol_3)
Sol_4 <- read.table(paste0(path, "01_Data/ForScmap/",
                           "E8.5_whole_embryo_Wnt1_counts.txt"), 
                    sep = " ")
Sol_4$gene <- rownames(Sol_4)
Sol_5 <- read.table(paste0(path, "01_Data/ForScmap/", 
                           "E9.5_anterior_Sox10_counts.txt"), 
                    sep = " ")
Sol_5$gene <- rownames(Sol_5)
Sol_6 <- read.table(paste0(path, "01_Data/ForScmap/",
                           "E9.5_posterior_Sox10_counts.txt"), 
                    sep = " ")
Sol_6$gene <- rownames(Sol_6)
Sol_7 <- read.table(paste0(path, "01_Data/ForScmap/",
                           "E9.5_trunk_Wnt1_counts.txt"), 
                    sep = " ")
Sol_7$gene <- rownames(Sol_7)

#Merge tables into one
Sol <- merge(Sol_1, Sol_2, by = "gene", all.x = TRUE, all.y = TRUE) %>%
  merge(Sol_3, by = "gene", all.x = TRUE, all.y = TRUE) %>%
  merge(Sol_4, by = "gene", all.x = TRUE, all.y = TRUE) %>%
  merge(Sol_5, by = "gene", all.x = TRUE, all.y = TRUE) %>%
  merge(Sol_6, by = "gene", all.x = TRUE, all.y = TRUE) %>%
  merge(Sol_7, by = "gene", all.x = TRUE, all.y = TRUE)

#Make nice
Sol[,-1] <- lapply(Sol[,-1], as.numeric)
Sol <- Sol[order(colnames(Sol))]

#Use gene names for rownames and remove gene column
rownames(Sol) <- Sol$gene
Sol <- dplyr::select(Sol, -c("gene"))

#Remove old tables
rm(Sol_1, Sol_2, Sol_3, Sol_4, Sol_5, Sol_6, Sol_7); gc()
#____________________________________________________________________________________
print("step 5")

## Annotations#______________________________________________________________________

#read data
Sol_ann <- read.table(paste0(path, "01_Data/ForScmap/", 
                             "emb_joint.txt"), 
                      sep = " ", comment.char = "", header = TRUE)
Sol_ann <- Sol_ann[order(Sol_ann$new.name), ]
rownames(Sol_ann) <- Sol_ann$new.name

#Subset only the cells that are in both the annotation and the count table
#Do this in both files.
Sol_ann <- subset(Sol_ann, new.name %in% colnames(Sol))
Sol <- Sol[, colnames(Sol) %in% Sol_ann$new.name]
#ncol(Sol) == nrow(Sol_ann)
gc()

#Prep annotation for scmap
Sol_ann <- Sol_ann["col"]
Sol_ann$col <- gsub("#49FF00FF", "Mesenchyme",     Sol_ann$col)
Sol_ann$col <- gsub("lightgrey", "NeuralTube",     Sol_ann$col)
Sol_ann$col <- gsub("#FF0000FF", "Autonomic",      Sol_ann$col)
Sol_ann$col <- gsub("#FFFF00FF", "MigratorySens",  Sol_ann$col)
Sol_ann$col <- gsub("#FF8000FF", "EarlyMigratory", Sol_ann$col)
Sol_ann$col <- gsub("#8000FFFF", "MigratoryAut",   Sol_ann$col)
Sol_ann$col <- gsub("#FF00DBFF", "Delaminatory",   Sol_ann$col)
Sol_ann$col <- gsub("#0092FFFF", "Sensory",        Sol_ann$col)
colnames(Sol_ann) <- "cell_type1"

#Sol_ann <- Sol_ann["col"]
#Sol_ann$col <- gsub("#49FF00FF", "Soldatov",     Sol_ann$col)
#Sol_ann$col <- gsub("lightgrey", "Soldatov",     Sol_ann$col)
#Sol_ann$col <- gsub("#FF0000FF", "Soldatov",      Sol_ann$col)
#Sol_ann$col <- gsub("#FFFF00FF", "Soldatov",  Sol_ann$col)
#Sol_ann$col <- gsub("#FF8000FF", "Soldatov", Sol_ann$col)
#Sol_ann$col <- gsub("#8000FFFF", "Soldatov",   Sol_ann$col)
#Sol_ann$col <- gsub("#FF00DBFF", "Soldatov",   Sol_ann$col)
#Sol_ann$col <- gsub("#0092FFFF", "Soldatov",        Sol_ann$col)
#colnames(Sol_ann) <- "cell_type1"


#Print celltype labels
unique(Sol_ann$cell_type1)
#____________________________________________________________________________________
print("step 6")

# Prep sce-objects
## Remove genes that are not in both data sets#______________________________________

Sol <- Sol[rownames(Sol) %in% rownames(Cao), ]
Cao <- Cao[rownames(Cao) %in% rownames(Sol), ]
gc()
#____________________________________________________________________________________
print("step 7")

## Make Cao sce-object#______________________________________________________________
#Inspect the reference data
#head(MCao_ann)
#head(Cao)[1:3, 1:3]

#Create SingleCellExperiment object
Cao_sce <- SingleCellExperiment(assays = list(normcounts = as.matrix(Cao)), 
                                colData = Cao_ann)
rm(Cao, Cao_ann); gc()
logcounts(Cao_sce) <- log2(normcounts(Cao_sce) + 1)
# use gene names as feature symbols
rowData(Cao_sce)$feature_symbol <- rownames(Cao_sce)
# remove features with duplicated names
Cao_sce <- Cao_sce[!duplicated(rownames(Cao_sce)), ]
#____________________________________________________________________________________
print("step 8")

## Make Soldatov sce-object#_________________________________________________________
#Inspect the query data
#Pmur[1:3, 1:3]
#head(Pmur_ann)

#Create SingleCellExperiment object
Sol_sce <- SingleCellExperiment(assays = list(normcounts = as.matrix(Sol)),
                                colData = Sol_ann)
rm(Sol, Sol_ann); gc()
#There are NAs in the matrix. replace those with 0s
assay(Sol_sce, "normcounts")[is.na(assay(Sol_sce, "normcounts"))] <- 0
#Logtransform
logcounts(Sol_sce) <- log2(normcounts(Sol_sce) + 1)
# use gene names as feature symbols
rowData(Sol_sce)$feature_symbol <- rownames(Sol_sce)
# remove features with duplicated names
Sol_sce <- Sol_sce[!duplicated(rownames(Sol_sce)), ]
gc()
#____________________________________________________________________________________
print("step 9")

# scmapCluster()
#Follow scmap vignet 
#https://bioconductor.org/packages/release/bioc/vignettes/scmap/inst/doc/scmap.html  

## Construct scmap cluster index#____________________________________________________

#Select the most informative features (genes)
Cao_sce <- selectFeatures(object = Cao_sce, 
                          n_features = 1500,
                          suppress_plot = TRUE) %>%
  #Calculate mean gene expression for every gene in each cluster
  indexCluster()
#table(rowData(Cao_sce)$scmap_features)

#head(metadata(Cao_sce)$scmap_cluster_index)
#heatmap(as.matrix(metadata(Cao_sce)$scmap_cluster_index))
#____________________________________________________________________________________
print("step 10")

## Project Soldatov on Cao#__________________________________________________________
#Calculate a good threshold
1/length(unique(Cao_sce$cell_type1)) #-> threshold should be >0.027

#Project Soldatov cells onto mouse cluster index
Sol_scClus_res <- scmapCluster(projection = Sol_sce,
                               index_list = list(Cao=metadata(Cao_sce)
                                                 $scmap_cluster_index),
                               threshold = 0.03)
#1/length(unique(Cao_sce$cell_type1)) #-> threshold should be >0.03

#Export results
saveRDS(Sol_scClus_res, paste0(path, 
                               "03_Results/scmapSolOnCao_FullData_dardel_output_", 
                               Sys.Date(),
                               ".Rds"))
#____________________________________________________________________________________
print("step 11")

#Plot the results#___________________________________________________________________

#head(Pmur_scClus_res$scmap_cluster_labs)
#head(Pmur_scClus_res$scmap_cluster_siml)
p1 <- SankeyForScmap(QuerryClusters = colData(Sol_sce)$cell_type1,
                     ReferenceCelltypes = Sol_scClus_res$scmap_cluster_labs[,'Cao'],
                     Title = "Soldatov on Cao", 
                     Subtitle = "scmapCluster()", 
                     Caption = paste0("#Soldatov cells =", dim(Sol_sce)[2], ". ",
                                      "#Cao cells =", dim(Cao_sce)[2], ". ",
                                      "n_feature = 1500",  ". ",
                                      "threshold = 0.03"))
p1
ggsave(plot = p1, filename = paste0(path,"03_Results/",
                                    "scmapCluster_SolOnCao_FullData_dardel_",
                                    Sys.Date(), ".svg"))
#____________________________________________________________________________________
print("step 12")

#Rename classifications and plot again#______________________________________________

Sol_scClus_res_renamed <- Sol_scClus_res
tempvec <- Sol_scClus_res_renamed$scmap_cluster_labs
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
Sol_scClus_res_renamed$scmap_cluster_labs <- tempvec
rm(tempvec)

p2 <- SankeyForScmap(QuerryClusters = colData(Sol_sce)$cell_type1,
                     ReferenceCelltypes = Sol_scClus_res_renamed$scmap_cluster_labs[,'Cao'],
                     Title = "Soldatov on Cao", 
                     Subtitle = "scmapCluster() on orig. mouse annot., renamed after proj.", 
                     Caption = paste0("#Soldatov cells =", dim(Sol_sce)[2], ". ",
                                      "#Cao cells =", dim(Cao_sce)[2], ". ",
                                      "n_feature = 1500",  ". ",
                                      "threshold = 0.03"))
p2
ggsave(plot = p2, filename = paste0(path,"03_Results/",
                                    "scmapCluster_SolOnCao_FullData_dardel_simp_",
                                    Sys.Date(), ".svg"))
#____________________________________________________________________________________
print("step 13")

#### Count cells_____________________________________________________________________
#How many were classified as NCC?
length(Sol_scClus_res_renamed$scmap_cluster_labs[Sol_scClus_res_renamed$scmap_cluster_labs ==
                                                    "NCC"])
length(Sol_scClus_res_renamed$scmap_cluster_labs[Sol_scClus_res_renamed$scmap_cluster_labs ==
                                                    "Not NCC"])
length(Sol_scClus_res_renamed$scmap_cluster_labs[Sol_scClus_res_renamed$scmap_cluster_labs ==
                                                    "Possible NCC"])
length(Sol_scClus_res_renamed$scmap_cluster_labs[Sol_scClus_res_renamed$scmap_cluster_labs ==
                                                    "unassigned"])
#____________________________________________________________________________________
print("step 14")

#Calculated runtime__________________________________________________________________
end_time <- Sys.time()
paste0("Elapsed time: ", end_time - start_time)
#____________________________________________________________________________________
