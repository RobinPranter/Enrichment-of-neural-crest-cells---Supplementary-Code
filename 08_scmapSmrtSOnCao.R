#Setup__________________________________________________________________________
options(future.globals.maxSize = 2 * 1024^4) # 2 TB 

rm(list = ls()); invisible(gc())
start_time <- Sys.time()

#Load libraries
library(Seurat)
library(magrittr)
library(SingleCellExperiment)
library(scmap)
library(ggplot2)
library(ggsankey)

path <- "~/WorkDirs/PmurNCCPilot_10X_rerun/"
#_______________________________________________________________________________

# Project data onto mouse atlas
##To predict NCCs in the data I will project it onto the mouse atlas from 
##Cao et al 2019 (DOI:10.1038/s41586-019-0969-x) using the program scmap 
##https://scmap.sanger.ac.uk/scmap/)


## Read gene orthologies________________________________________________________
#Read file
OrthoTab <-read.table(paste0(path, 
                             "01_Data/BM_GCRm39vsPodMur1.0_ManDL240213.txt"), 
                      sep = ",", header = TRUE, na.strings = "") %>% na.omit() 
#This file was downloaded 2024-02-13 from Ensembl using BioMart with the 
#following settings:
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
#_______________________________________________________________________________
print("Step 1")

## Read mouse data
### Read count table____________________________________________________________
Mmus <- readRDS(paste0(path, "01_Data/CaoData/gene_count_cleaned.RDS")) %>% 
  as.matrix() %>% as.data.frame() #Transform sparse matrix to data frame
gc()

#Remove version number from Ensembl gene names
rownames(Mmus) <- gsub(pattern = "\\.\\d+$", 
                       replacement = "",
                       x = rownames(Mmus))
#_______________________________________________________________________________
print("Step 2")

#### Rename genes_______________________________________________________________
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
#_______________________________________________________________________________
print("Step 3")

### Read annotation file________________________________________________________
#Read data
Mmus_ann <- read.csv(paste0(path, "01_Data/CaoData/cell_annotate.csv")
)[ , c("sample", "Main_cell_type")]
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
#_______________________________________________________________________________
print("Step 4")

## Read and prep lizard data____________________________________________________
#Read the data
filename <- list.files(paste0(path, "02_Analysis/"), 
                       pattern = "02_FilterSmrtS_output_")
if (length(filename) == 1) {
  #Count table
  filt.CT <- readRDS(file = paste0(path, "02_Analysis/", filename)) %>%
    JoinLayers(assay = "RNA") %>%
    LayerData(assay = "RNA", layer = "counts") %>% 
    as.data.frame()
  #Annotations
  Pmur_ann <- readRDS(file = paste0(path, "02_Analysis/", filename)
  )[["sampleID"]] %>%as.data.frame()
} else {
  print("Choose file to read")
  print(filename)
}
gc()

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
gc()
#_______________________________________________________________________________
print("Step 8")

## scmapCluster()
#Follow scmap vignet 
#https://bioconductor.org/packages/release/bioc/vignettes/scmap/inst/doc/scmap.html 

### Remove genes that are not in both data sets_________________________________
filt.CT <- filt.CT[rownames(filt.CT) %in% rownames(Mmus), ]
Mmus <- Mmus[rownames(Mmus) %in% rownames(filt.CT), ]
gc()
#_______________________________________________________________________________
print("Step 9")

### Create lizard sce___________________________________________________________
#Creating the single cell experiment object for the lizard data because that is
#the format that scmap takes.

#Inspect the query data
#filt.CT[1:3, 1:3]
#head(Pmur_ann)

#Create SingleCellExperiment object
Pmur_sce <- SingleCellExperiment(assays = list(normcounts = as.matrix(filt.CT)),
                                 colData = Pmur_ann)
rm(Pmur_ann, filt.CT); gc()
logcounts(Pmur_sce) <- log2(normcounts(Pmur_sce) + 1)
# use gene names as feature symbols
rowData(Pmur_sce)$feature_symbol <- rownames(Pmur_sce)
# remove features with duplicated names
Pmur_sce <- Pmur_sce[!duplicated(rownames(Pmur_sce)), ]

gc()
#_______________________________________________________________________________
print("Step 10")

### Create mouse sce____________________________________________________________
#Creating the single cell experiment object for the mouse data because that is
#the format that scmap takes.
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
#_______________________________________________________________________________
print("Step 11")

### Construct scmap cluster index_______________________________________________
#Select the most informative features (genes)
Mmus_sce <- selectFeatures(object = Mmus_sce, 
                           n_features = 1500,
                           suppress_plot = TRUE) %>% 
  #Calculate mean gene expression for every gene in each cluster
  indexCluster()

#head(metadata(Mmus_sce)$scmap_cluster_index)
#heatmap(as.matrix(metadata(Mmus_sce)$scmap_cluster_index))
#_______________________________________________________________________________
print("Step 12")
#### Project lizard on mouse____________________________________________________
#Calculate a good threshold
1/length(unique(Mmus_sce$cell_type1)) #-> threshold should be >0.027

#Project lizard cells onto mouse cluster index
Pmur_scClus_res <- scmapCluster(projection = Pmur_sce,
                                index_list = list(Cao=metadata(Mmus_sce)$
                                                    scmap_cluster_index),
                                threshold = 0.03)
#head(Pmur_scClus_res$scmap_cluster_labs)
#head(Pmur_scClus_res$scmap_cluster_siml)
#_______________________________________________________________________________
print("Step 13")

#Save result data_______________________________________________________________
saveRDS(list(Pmur_sce, Pmur_scClus_res), 
        paste0(path, "02_Analysis/06_scmapSmrtSOnCao_output_",
               Sys.Date(), ".pdf"))
#_______________________________________________________________________________
print("Step 14")

#### Inspect the results________________________________________________________
#Define function for plotting
SankeyForScmap <- function(QuerryClusters, 
                           ReferenceCelltypes,
                           Title = "HereGoesATitle",
                           Subtitle = "HereGoesASubtitle",
                           Caption = "HereGoesACaption"){
  #Setup dataframe
  df <- data.frame(cbind(QuerryClusters, ReferenceCelltypes))
  names(df) <- c('QuerryClusters', 'ReferenceCellTypes')
  df <- df %>% ggsankey::make_long(QuerryClusters, ReferenceCellTypes)
  
  #Construct plot
  pl <- ggplot(df, aes(x = x,                        
                       next_x = next_x,                                     
                       node = node,
                       next_node = next_node,        
                       fill = factor(node),
                       label = node)) 
  pl <- pl + geom_sankey(flow.alpha = 0.5,     
                         node.color = "black",
                         show.legend = TRUE)
  
  pl <- pl + geom_sankey_text(color = "black") 
  pl <- pl + theme_bw()
  pl <- pl + theme(legend.position = 'none')
  pl <- pl + theme(axis.title = element_blank(),
                   axis.text.y = element_blank(),
                   axis.ticks = element_blank(),
                   panel.grid = element_blank())
  
  pl <- pl + scale_fill_viridis_d(option = "inferno")
  pl <- pl + labs(title = Title)
  pl <- pl + labs(subtitle = Subtitle)
  pl <- pl + labs(caption = Caption)
  pl <- pl + labs(fill = 'Nodes')
  
  pl
}

p1 <- SankeyForScmap(QuerryClusters = colData(Pmur_sce)$sampleID,
                     ReferenceCelltypes = Pmur_scClus_res$scmap_cluster_labs[,'Cao'],
                     Title = "scmap SmrtS on Cao",
                     Subtitle = "scmapCluster()",
                     Caption = paste0("#Lizard cells=", dim(Pmur_sce)[2], ". ",
                                      "#Mouse cells=", dim(Mmus_sce)[2], ". ",
                                      "n_feature=1500. ",
                                      "threshold=0.03."))

ggsave(plot = p1, filename = paste0(path,"03_Results/", 
                                    "06_scmapSmrtSOnCao_",
                                    Sys.Date(), ".pdf"))
#_______________________________________________________________________________
print("Step 15")

#### Rename classifications and plot again______________________________________
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
                     Title = "scmap SmrtS on Cao",
                     Subtitle = "scmapCluster() renamed after proj",
                     Caption = paste0("#Lizard cells=", dim(Pmur_sce)[2], ". ",
                                      "#Mouse cells=", dim(Mmus_sce)[2], ". ",
                                      "n_feature=1500. ",
                                      "threshold=0.03."))
ggsave(plot = p2, filename = paste0(path,"03_Results/",
                                    "06_scmapSmrtSOnCao_simp_",
                                    Sys.Date(), ".pdf"))
#_______________________________________________________________________________
print("Step 16")
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
print("Step 17")

#Calculated runtime_____________________________________________________________
end_time <- Sys.time()
paste0("Elapsed time: ", end_time - start_time)
#_______________________________________________________________________________


























