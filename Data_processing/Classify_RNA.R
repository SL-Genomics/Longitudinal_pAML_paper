#script to classify clusters within scRNA data for each patient sample
#This script is meant to illustrate the broad workflow, in principle any number of samples can be combined in Preprocess_RNA
#This workflow needs to be adjusted appropriately to accomodate multiple samples or multiple patients
#In the study typically diagnosis remission and relapse samples were combined for each analysis

path_annotated_scRNA <- #paths of output locations for each patient from Annotate_RNA.R
path_scRNA <- #path of output location of processed scRNA data for each patient

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

# run Annotate_RNA.R first to acquire this
scRNA_data<- readRDS(path_annotated_scRNA_out)

scRNA_data <- FindNeighbors(RNA_cells_anno, graph.name = "NN")
scRNA_data <- FindClusters(RNA_cells_anno, graph.name = "NN")

#classification based on mutations

cell_ID_clust<-list()
for (clust in str_sort(unique(scRNA_data$seurat_clusters),numeric=T))
{
    cell_ID_clust[[clust]]<-table(scRNA_data[scRNA_data$seurat_clusters == clust,]$mal_binary)
}
all_IDs<-do.call(rbind,cell_ID_clust)

if (length(which(scRNA_data$mal_binary %in% "Malignant")) >= 1)
{
fishers<-list()
#do fisher tests
for (clust in str_sort(unique(df_BIO_ID$Clusters),numeric=T))
{
p.table <- data.frame(
  "Malignant" = c(all_IDs[clust,"Malignant"],sum(all_IDs[rownames(all_IDs)[rownames(all_IDs) %ni% clust],"Malignant"]) ),
  "Not_mal" = c(all_IDs[clust,"No malignant allele"], sum(all_IDs[rownames(all_IDs)[rownames(all_IDs) %ni% clust],"No malignant allele"])),
  row.names = c("Cluster", "Not_cluster"),
  stringsAsFactors = FALSE
)
fishers[[clust]]<-p.adjust(fisher.test(p.table,alternative = "greater")$p.value,method="BH",n=length(unique(scRNA_data$seurat_clusters)))
}
}

#manually assign clusters based on annotations and p-values
# i.e. Malignant<-c(0,1,2,3,7,12,20)

scRNA_data$Malignant_cluster<- "Normal"
scRNA_data$Malignant_cluster[RNA_cells_anno$seurat_clusters %in% Malignant]<-"Malignant"

saveRDS(scRNA_data,path_scRNA)
