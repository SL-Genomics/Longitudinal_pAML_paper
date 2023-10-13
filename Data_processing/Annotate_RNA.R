#script to add annotations for cell type classification, identified mutations and doublets
#This script is meant to illustrate the broad workflow, in principle any number of samples can be combined in Preprocess_RNA
#This workflow needs to be adjusted appropriately to accomodate multiple samples or multiple patients
#In the study typically diagnosis remission and relapse samples were combined for each analysis

path_reference_out<- #path to reference data generated using Process_reference.R
path_processed_scRNA_out<- #paths of output locations for each patient from Preprocess_RNA.R
CB_sniffer_out<- #path to output (tsv) of the CB_snifffer tool 
path_annotated_scRNA <- #paths of output locations for each patient
path_classification_stats<-#path to output predictions for each celltype

#load libraries
require(dplyr)
require(Seurat)
require(patchwork)
require(hdf5r)
require(DropletUtils)
require(scPred)
require(magrittr)
require(SeuratData)
require(ScDblFinder)

# to make the script reproducible
set.seed(42)

# run Preprocess_RNA.R first to acquire this
scRNA_data<- readRDS(path_processed_scRNA_out)

# run Process_reference.R first to acquire this
reference<- readRDS(path_reference_out)

#run scPred on dataset and output a clustering based on the model
sc_pred_query <- scRNA_data
sc_pred_query <- scPredict(sc_pred_query, reference)
sc_pred_query <- RunUMAP(sc_pred_query, reduction = "scpred", dims = 1:30)
write.table(sc_pred_query@meta.data,file = path_classification_stats,append = F,sep = "\t",row.names=F)

# apply cell predictions to original clustering with and without rejection (low scoring classifications)
scRNA_data$pred<-sc_pred_query$scpred_prediction
scRNA_data$pred_no_reject<-sc_pred_query$scpred_no_rejection

# add malignant classification to cells based on cb sniffer
bc_counts<-read.delim(CB_sniffer_out)

barcode_alt<-list()
for (cb in rownames(scRNA_data))
{	
barcode_alt<-sapply(rownames(scRNA_data),function(cb){mean(bc_counts[bc_counts$barcode %in% cb,]$alt_count))
}
barcode_alt[is.na(barcode_alt)]<--1

scRNA_data$CB_score<-barcode_mat[rownames(scRNA_data),]
scRNA_data$mal_binary<-scRNA_data$CB_score 

scRNA_data$mal_binary[scRNA_data$mal_binary > 0] <- "Malignant_allele"
scRNA_data$mal_binary[scRNA_data$mal_binary == -1] <- "Unknown"
scRNA_data$mal_binary[scRNA_data$mal_binary == 0] <- "No_malignant_allele"

# add doublet scores to each cell using scDBLfinder
scRNA_data_sce<-as.SingleCellExperiment(scRNA_data)
scRNA_data_dbl<-scDblFinder(scRNA_data_sce,samples=scRNA_data_sce$orig.ident)
scRNA_data<-as.Seurat(scRNA_data_dbl)
scRNA_data$db_density<-log2(scDblFinder::computeDoubletDensity(scRNA_data_sce)+1)

saveRDS(scRNA_data,path_annotated_scRNA)