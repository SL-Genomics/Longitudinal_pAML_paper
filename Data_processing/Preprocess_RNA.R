#script to preprocess counts for each patient sample
#This script is meant to illustrate the broad workflow, in principle any number of samples can be combined
#This workflow needs to be adjusted appropriately to accomodate multiple samples or multiple patients
#In the study typically diagnosis remission and relapse samples were combined for each analysis

paths_cellranger_out<- #paths to count matrices out of cellranger for each patient
path_processed_scRNA_out<- # paths of output locations for each sample (i.e. diagnosis, relapse, remission)
barcodes_out<- #path to write a table of filtered barcodes (necessary for CB_sniffer)

#load libraries
require(dplyr)
require(Seurat)
require(patchwork)
require(hdf5r)
require(DropletUtils)
require(scPred)
require(magrittr)
require(SeuratData)

# to make the script reproducible
set.seed(42)

# combine samples and add barcodes to the rownames for easy filtering
merge_mat <- read10xCounts(samples= paths_cellranger_out, col.names=TRUE, type="sparse")
colnames(merge_mat) <- paste0(merge_mat$Sample,"_",merge_mat$Barcode)
rownames(merge_mat) <- rowData(merge_mat)$Symbol
merge_mat.s <- as(counts(merge_mat), "dgCMatrix")

# Seurat normalization and filtering (> 1000 UMI, > 500 genes, 15% MT reads)
scRNA_data <- Seurat::CreateSeuratObject(counts = merge_mat.s, project = "pAML", metadata= merge_mat$Sample)
scRNA_data<-subset(scRNA_data,subset= nCount_RNA > 1000)  
scRNA_data<-subset(scRNA_data,subset= nFeature_RNA > 500)
scRNA_data[["percent.mt"]] <- PercentageFeatureSet(scRNA_data, pattern = "^MT-")
scRNA_data <- subset(scRNA_data, subset = percent.mt < 15)
scRNA_data <- scRNA_data %>% 
  NormalizeData() %>% 
  FindVariableFeatures() %>% 
  ScaleData() %>% 
  RunPCA() %>% 
  RunUMAP(dims = 1:30)
  
write.table(rownames(scRNA_data@meta.data), barcodes_out)  
saveRDS(scRNA_data,path_processed_scRNA)








