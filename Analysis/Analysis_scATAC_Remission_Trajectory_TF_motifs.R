# script to analyze generate a reference of remission samples to infer trajectories
# requires processed scATAC data (Classify_ATAC.R) to work 

path_scATAC_remission <- #path of output location of combined processed scATAC data for all patients (Classify_ATAC.R), this could be a normal bone marrow reference as well
outputdir<- #path to output the peak set used for lineage priming prediction

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
require(GenomicRanges)

# to make the script reproducible
set.seed(42)

scATAC_data<-loadArchRProject(path_scATAC_remission)

# first test where celltypes are located in the clustering
p_remission<- plotEmbedding(ArchRProj = scATAC_data, colorBy = "cellColData", name = "scpred_no_rejection", embedding = "UMAP",pal = color_palette_celltypes,randomize = F,rastr=F)
p_cluster<- plotEmbedding(ArchRProj = scATAC_data, colorBy = "cellColData", name = "Clusters", embedding = "UMAP",randomize = F,rastr=F)
plotPDF(p_remission,p_cluster name = "Celltypes_UMAP.pdf", ArchRProj = scATAC_data, addDOC = FALSE, width = 5, height = 5)	

#check for identity of clusters (identified in preprocess_atac) to identify logical lineages, clusters work better than celltype predictions unless the clustering is really clean
#i.e. plot celltype proportions per cluster and plot marker gene expression/ TF motif ernichment
#i.e. myeloidclusters<- c(cluster enriched in early cells(HSC/progenitor), cluster(s) in between(GMP), differentiated cluster(Monocyte))

scATAC_data<- addTrajectory(
	ArchRProj = scATAC_data, 
	name = Trajectory_name, 
	groupBy = "Clusters",
	trajectory = clusters, 
	embedding = "UMAP", 
	force = TRUE,
)

p_trajectory <- plotTrajectory(scATAC_data, trajectory = Trajectory_name, colorBy = "cellColData", name = Trajectory_name)
plotPDF(p_remission,p_trajectory, name = "Celltypes_Trajectory.pdf", ArchRProj = scATAC_data, addDOC = FALSE, width = 5, height = 5)

#try different combinations to identify a sensible lineage

trajectory_motif  <- getTrajectory(ArchRProj = scATAC_data, name = Trajectory_name, useMatrix = "MotifMatrix", log2Norm = FALSE)
trajectory_expression <- getTrajectory(ArchRProj = scATAC_data, name = Trajectory_name, useMatrix = "GeneIntegrationMatrix", log2Norm = TRUE)

#use correlated motif and gene expression trajectories only
correlated_trajectory <- correlateTrajectories(trajectory_expression, trajectory_motif)

all_names<- unique(correlated_trajectory$correlatedMappings$matchname1)
matchlist<-correlated_trajectory[[2]][correlated_trajectory[[2]]$matchname1 %in% all_names,]

matchlist<-matchlist[1:nrow(matchlist) %ni% grep("deviations",matchlist$name2),]
matchlist<-matchlist[1:nrow(matchlist) %ni% grep("-",matchlist$name1),]

trajectory_motif_subset <- trajectory_motif[matchlist$name2, ]
trajectory_expression_subset <- trajectory_expression[matchlist$name1, ]

rownames(trajectory_motif_subset)<-gsub('_[[:digit:]]+',"",as.vector(sapply(trajectory_motif_subset,function(x){unlist(strsplit(rownames(x),"z:",fixed=T))[2]})))
rownames(trajectory_expression_subset)<-as.vector(sapply(trajectory_expression_subset,function(x){unlist(strsplit(rownames(x),":",fixed=T))[2]}))

combinedMat <- plotTrajectoryHeatmap(trajectory_motif_subset, returnMat = TRUE, varCutOff = 0)
ordered_rows<-rownames(combinedMat)

ht1 <- plotTrajectoryHeatmap(trajMM2, pal = paletteContinuous(set = "solarExtra"), varCutOff = 0,rowOrder= ordered_rows)
ht2 <- plotTrajectoryHeatmap(trajGIM2, pal = paletteContinuous(set = "blueYellow"), varCutOff = 0, rowOrder= ordered_rows)

plotPDF(ht1 + ht2, name = "Trajectory_Heatmaps.pdf", ArchRProj = scATAC_data, addDOC = FALSE, width = 10, height = 8)




