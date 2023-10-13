# script to analyze generate a reference of remission samples to infer trajectories and output regions used to infer priming
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

trajectory_Peaks  <- getTrajectory(ArchRProj = scATAC_data, name = Trajectory_name, useMatrix = "PeakMatrix", log2Norm = TRUE)
trajectory_peak_matrix <- plotTrajectoryHeatmap(trajectory_Peaks, pal = paletteContinuous(set = "solarExtra"),returnMatrix = T)

#add cluster annotations to annotate lineages across pseudotime, ArchR divides the trajectory in 100 bins

df<-cbind(getCellColData(scATAC_data)[[Trajectory_name]],getCellColData(scATAC_data)$Clusters)   
df<-df[!is.na(df[,1]),]
anno<-df[order(df[,1]),2]

# get the majority of names (could be clusters or something else including celltypes) present in each bin
anno_squeeze<-vector()
for(bin in 0:99)
{
    anno_squeeze[bin+1]<-names(which.max(table(anno[((bin*floor(length(anno)/100))+1):((bin+1)*floor(length(anno)/100))])))
}

#Identify specific peaks in each trajectory used for lineage priming scores
#use the most primitive cluster (stem_cluster) to compare which peaks are upregulated

mean_score_peaks<-mean(trajectory_peak_matrix)
Trajectory_peaks<-trajectory_peak_matrix[apply(trajectory_peak_matrix[,which(anno_squeeze != stem_cluster)],1,median) >mean_score_peaks,]

chroms<-sapply(rownames(Trajectory_peaks),function(x){unlist(strsplit(x,split=":",fixed=T))[1]})
locs<-sapply(rownames(Trajectory_peaks),function(x){unlist(strsplit(x,split=":",fixed=T))[2]})
starts<-sapply(locs,function(x){unlist(strsplit(x,split="_",fixed=T))[1]})
ends<-sapply(locs,function(x){unlist(strsplit(x,split="_",fixed=T))[2]})
variance<-apply(Trajectory_peaks,1,var)
ranged_object<-GRanges(seqnames=chroms,ranges=IRanges(start = as.numeric(starts),end=as.numeric(ends)),strand = "*",score=variance)
seqlevels(ranged_object)<-str_sort(x = seqlevels(ranged_object),numeric=T)
Trajectory_ranges<-sort(ranged_object)

#bgd_peaks are peaks found in other trajectories i.e. lymphoid peaks, this is to make sure peaks are unique to avoid lineage ambiguity
top_peaks<-5000
Trajectory_ranges_select<-Trajectory_ranges[1:length(Trajectory_ranges) %ni% queryHits(findOverlaps(Trajectory_ranges,bgd_peaks,minoverlap = 1)),]
Trajectory_ranges_select<-Trajectory_ranges_select[order(Trajectory_ranges_select$score)<=top_peaks,]

write.table(Trajectory_ranges_select,outputdir,sep="\t",quote=F,row.names=F)



