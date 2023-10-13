# script to filter scATAC data and run basic peak calling
#This script is meant to illustrate the broad workflow, in principle any number of samples can be combined
#This workflow needs to be adjusted appropriately to accomodate multiple samples or multiple patients
#In the study typically diagnosis remission and relapse samples were combined for each analysis

paths_cellranger_atac_out <- #paths to fragment files from cellranger_atac for each sample, as a named vector
output_dir <- #path of output location of processed scATAC data for each patient
MACS2_location <- #path to MACS2 (https://pypi.org/project/MACS2/)
output_dir_barcodes<- #path to save barcodes 

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

#set genome to hg38
addArchRGenome("hg38")

dir.create(output_dir)
setwd(output_dir)

#read files and turn into arrowsfiles
#also generates QC 	
#inputFiles can be any number of files (i.e. diagnosis remissiona and relapse or remissions across patients)
ArrowFiles <- createArrowFiles(inputFiles = paths_cellranger_atac_out, sampleNames = names(paths_cellranger_atac_out), 
minTSS = 4, minFrags = 2000, addTileMat = T, addGeneScoreMat = T)

# load in the data as arrowfiles	
scATAC_data <- ArchRProject(ArrowFiles = ArrowFiles, outputDirectory = output_dir)

# filter for % reads in TSS, works better for me than an arbitrary TSS enrichment value	
metadata <-getCellColData(scATAC_data)
TSS_ratio<-metadata$ReadsInTSS/metadata$Total_read_number
metadata<-metadata[TSS_ratio >= 0.1,]
	
# remove cells that do not pass the threshold				
scATAC_data <-subsetCells(scATAC_data,rownames(metadata))

# Add scores for doublet calling and remove putative doublets
doubletscores <- addDoubletScores(input = scATAC_data, k = 10, knnMethod = "UMAP",LSIMethod = 1)
scATAC_data <- filterDoublets(ArchRProj = doubScores)	
	
#perform clustering
scATAC_data <- addIterativeLSI(ArchRProj = scATAC_data, useMatrix = "TileMatrix", name = "IterativeLSI",force=T)
scATAC_data <- addClusters(input = scATAC_data, reducedDims = "IterativeLSI",force=T)
scATAC_data <- addUMAP(ArchRProj = scATAC_data, reducedDims = "IterativeLSI",force=T)
scATAC_data <- addTSNE(ArchRProj = scATAC_data, reducedDims = "IterativeLSI", name = "TSNELSI", perplexity = 30,force=T)

# write a table with all included barcodes at this point, used for filtering in CB sniffer
barcodes<-rownames(getCellColData(scATAC_data))
write.table(barcodes,output_dir)

#call and add peaks and infer motif enrichment
scATAC_data <- addImputeWeights(scATAC_data)
scATAC_data <- addGroupCoverages(ArchRProj = scATAC_data, groupBy = "Clusters", force= T)
scATAC_data <- addReproduciblePeakSet(ArchRProj = scATAC_data, groupBy = "Clusters", pathToMacs2 = MACS2_location, 
minCells = 50,method = "p",cutOff = 1e-5,peakMethod = "Macs2")	
scATAC_data <- addPeakMatrix(scATAC_data)
scATAC_data <- addMotifAnnotations(ArchRProj = scATAC_data, motifSet = "cisbp", name = "Motif",force=T)
scATAC_data <- addBgdPeaks(scATAC_data)
scATAC_data <- addDeviationsMatrix(ArchRProj = scATAC_data, peakAnnotation = "Motif", force = TRUE) 
scATAC_data <- addArchRAnnotations(ArchRProj = scATAC_data, collection = "ATAC", force = TRUE)

#save RDS file in outputdirectory	
saveArchRProject(ArchRProj = scATAC_data, outputDirectory = output_dir, load = FALSE)