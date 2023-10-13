#script to classify clusters within scATAC data for each patient sample
#This script is meant to illustrate the broad workflow, in principle any number of samples can be combined in Preprocess_ATAC 
#This workflow needs to be adjusted appropriately to accomodate multiple samples or multiple patients
#In the study typically diagnosis remission and relapse samples were combined for each analysis

output_dir <- #paths of output locations for each patient from Annotate_ATAC.R

#load libraries
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

# run Annotate_RNA.R first to acquire this
scATAC_data <-loadArchRProject(output_dir)

#classification based on mutations

metadata<-getCellColData(scATAC_data)

cell_ID_clust<-list()
for (clust in str_sort(unique(metadata$Clusters),numeric=T))
{
    cell_ID_clust[[clust]]<-table(metadata[metadata$Clusters == clust,]$mal_binary)
}
all_IDs<-do.call(rbind,cell_ID_clust)

if (length(which(metadata$mal_binary %in% "Malignant")) >= 1)
{
fishers<-list()
#do fisher tests
for (clust in str_sort(unique(metadata$Clusters),numeric=T))
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

metadata$Malignant_cluster<- "Normal"
metadata$Malignant_cluster[metadata$Clusters %in% Malignant]<-"Malignant"

scATAC_data<-addCellColData(ArchRProj=scATAC_data,data=metadata$Malignant_cluster,name = "Malignant_cluster",cells = scATAC_data$cellNames,force=T)

#save RDS file in outputdirectory	
saveArchRProject(ArchRProj = scATAC_data, outputDirectory = output_dir, load = FALSE)
