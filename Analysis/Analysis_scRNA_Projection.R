# script to project cells upon a reference and to identify differences in abundance 
# requires process_reference.R and processed scRNA data (Classify_RNA.R) to work 

path_scRNA <- #path of output location of processed scRNA data for each patient (Classify_RNA.R)
path_reference_out<- #path to processed reference (Process_Reference.R)

#load libraries
require(dplyr)
require(Seurat)
require(patchwork)
require(hdf5r)
require(DropletUtils)
require(scPred)
require(magrittr)
require(SeuratData)
require(symphony)
require(miloR)
require(SingleCellExperiment)
require(scater)
require(scran)

#requires functions from symphony to be imported (https://github.com/immunogenomics/symphony/blob/main/vignettes/utils_seurat.R)
source('utils_Seurat.R')

#make sure this is the same as process reference for references to be comparable
set.seed(42)

#load in RNA data, make sure to subset data for the population of interest (i.e. malignant vs non-malignant)
scRNA_data<- readRDS(path_scRNA)

#reprocess the reference to include harmony before running UMAP and save the model 
harm_ref<-ref %>% 
NormalizeData() %>% 
FindVariableFeatures() %>% 
ScaleData() %>% 
RunPCA() %>% 
RunHarmony.Seurat('orig.ident') %>% 
FindNeighbors(dims = 1:20, reduction = 'harmony') %>% 
FindClusters(resolution = 0.5)
   
harm_ref<- RunUMAP(Embeddings(harm_ref, 'harmony')[, 1:20], assay='RNA', verbose=FALSE, umap.method='uwot', return.model=TRUE)
harm_ref <- buildReferenceFromSeurat(harm_ref, verbose = TRUE, save_umap = TRUE, save_uwot_path = 'cache_symphony.uwot')

#make sure the reference looks as expected
DimPlot(harm_ref, reduction = 'umap', shuffle = TRUE)

#map cells on the reference
Symphony_data <- mapQuery(scRNA_data@assays$RNA@counts, scRNA_data@meta.data, harm_ref,vars = 'orig.ident', return_type = 'Seurat')

#save coordinates for plotting
coords <- Embeddings(Symphony_data[["umap"]])
coords_ref <- Embeddings(harm_ref[["umap"]])

#example of plotting function
ggplot(data=as.data.frame(coords_ref),aes(x=UMAP_1,y=UMAP_2)) +
geom_point(color="#EEEEEE",shape=16)+
geom_hex(data = as.data.frame(coords),aes(x=UMAP_1,y=UMAP_2,alpha=(..count..)),color="#99999911",bins=20) +
scale_fill_gradient(trans="sqrt",low = "grey",high="#111111")+ scale_alpha_continuous(trans="sqrt",range = c(0,1))+
theme_minimal()

#compare abundance

Symphony_data_sce <- as.SingleCellExperiment(Symphony_data)
Symphony_data_milo <- Milo(Symphony_data_sce)

Milo_object <- buildGraph(Symphony_data_milo, k = 100, d = 30)
Milo_object <- makeNhoods(Milo_object, prop = 0.1, k = 100, d=30, refined = TRUE)

#check the NHood distribution to obtain valuable sizes
plotNhoodSizeHist(Milo_object)

Milo_object <- countCells(Milo_object, meta.data = data.frame(colData(Milo_object)), samples="orig.ident")

#Biopsy is a variable containing cell biopsy origin (i.e. Diagnosis and Relapse)
Milo_object_design <- data.frame(colData(Milo_object))[,c(orig.ident,Biopsy)]
Milo_object_design <- distinct(Milo_object_design)
rownames(Milo_object_design) <- Milo_object_design$orig.ident

## Reorder rownames to match columns 
Milo_object_design <- Milo_object_design[colnames(nhoodCounts(Milo_object_design)), , drop=FALSE]

Milo_object <- calcNhoodDistance(Milo_object, d=30)
rownames(Milo_object_design) <- Milo_object_design$orig.ident
DA_results <- testNhoods(Milo_object, design = ~ Biopsy, design.df = Milo_object_design)
#Biopsy is a variable containing cell classification (i.e. HSC and monocyte)
DA_results <- annotateNhoods(Milo_object, DA_results, coldata_col = Cell_type)
Milo_object_design <- buildNhoodGraph(Milo_object,overlap = 1000)




