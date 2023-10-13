# script to link scRNA and scATAC data, perform peak to gene linkage and to annotate cells
#This script is meant to illustrate the broad workflow, in principle any number of samples can be combined in Preprocess_ATAC
#This workflow needs to be adjusted appropriately to accomodate multiple samples or multiple patients
#In the study typically diagnosis remission and relapse samples were combined for each analysis

output_dir <- #path of output location of processed scATAC data for each patient after running Preprocess_ATAC.R
path_annotated_scRNA <- #paths of output locations for each patient from Annotate_RNA.R
CB_sniffer_out_ATAC<- #path to output (tsv) of the CB_sniffer tool performed on scATAC data

# load libraries
require(ArchR)
require(Seurat)
require(TxDb.Hsapiens.UCSC.hg38.knownGene)
require(BSgenome.Hsapiens.UCSC.hg38)
require(magrittr)
require(dplyr)
require(patchwork)
require(stringr)
require(hdf5r)
require(DropletUtils)

# to make the script reproducible
set.seed(42)

setwd(output_dir)

# run Annotate_RNA.R first to acquire this
scRNA_data <- readRDS(path_processed_scRNA_out)

#run this script only after finishing Preprocess_ATAC.R
scATAC_data <-loadArchRProject(output_dir)

#link scRNA data to scATAC data
scATAC_data <- addGeneIntegrationMatrix(ArchRProj = proj, useMatrix = "GeneScoreMatrix", matrixName = "GeneIntegrationMatrix", 
reducedDims = "IterativeLSI",seRNA = scRNA_data, addToArrow = T, groupRNA = "orig.ident", force = T)

#transfer the celltype labels from the linked cells
celltypes<-RNA@meta.data
metadata<-getCellColData(scATAC_data)
metadata$scpred_prediction<-celltypes[match(metadata$predictedCell,rownames(celltypes)),]$scpred_prediction
metadata$scpred_no_rejection<-celltypes[match(metadata$predictedCell,rownames(celltypes)),]$scpred_no_rejection

# add malignant classification to cells based on cb sniffer
bc_counts<-read.delim(CB_sniffer_out_ATAC)

barcode_alt<-list()
for (cb in rownames(metadata))
{	
barcode_alt<-sapply(rownames(metadata),function(cb){mean(bc_counts[bc_counts$barcode %in% cb,]$alt_count))
}
barcode_alt[is.na(barcode_alt)]<--1

metadata$CB_score<-barcode_mat[rownames(metadata),]
metadata$mal_binary<-metadata$CB_score 

metadata$mal_binary[metadata$mal_binary > 0] <- "Malignant_allele"
metadata$mal_binary[metadata$mal_binary == -1] <- "Unknown"
metadata$mal_binary[metadata$mal_binary == 0] <- "No_malignant_allele"

# Add annotations to cells
for (idx in c("CB_score","mal_binary","scpred_prediction","scpred_no_rejection"))
{
scATAC_data<-addCellColData(ArchRProj=scATAC_data,data=metadata[,idx],name = idx,cells = scATAC_data$cellNames,force=T)
}

#save RDS file in outputdirectory	
saveArchRProject(ArchRProj = scATAC_data, outputDirectory = output_dir, load = FALSE)
