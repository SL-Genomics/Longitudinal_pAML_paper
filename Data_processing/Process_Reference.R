#script to train a model for cell type identification
#run once, run manually, will be different for every reference

path_ref_in<- #path to reference data containing annotations
path_reference_out<-#output path

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

#load in the reference
ref <-readRDS(path_ref_in)

#process standard reference and turn into seurat
refs <- as(bm, "SingleCellExperiment")
ref <- CreateSeuratObject(counts = counts(bms), meta.data = as.data.frame(colData(bms)))

#filter reads in the same way as tumor data
ref<-subset(ref,subset= nCount_RNA > 1000)  
dat<-subset(ref,subset= nFeature_RNA > 500)
ref[["percent.mt"]] <- PercentageFeatureSet(ref, pattern = "^MT-")
ref <- subset(ref, subset = percent.mt < 15)

#normalize and run UMAP
ref <- ref %>% 
NormalizeData() %>% 
FindVariableFeatures() %>% 
ScaleData() %>% 
RunPCA() %>% 
RunUMAP(dims = 1:30) 

#check whether the classifications match with what was published
DimPlot(object = reference,group.by = "BioClassification",shuffle = T)

#Train model
reference <- getFeatureSpace(reference, "BioClassification")
reference <- trainModel(reference)
get_scpred(reference)

#check the model for false positive and false negative classifications using the input
plot_probabilities(reference)

#save for easy loading
saveRDS(reference,path_reference_out)







