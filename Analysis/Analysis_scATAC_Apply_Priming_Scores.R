# script to analyze priming scores within HSPCs, requires Analysis_scATAC_Remission_Trajectory.R to run
# requires processed scATAC data (Classify_ATAC.R) and a processed trajectory on normal cells (Analysis_scATAC_Remission_Trajectory_priming.R) to work 

path_scATAC <- #path of output location of processed scATAC data for a comparison (i.e. HSPCs) (Classify_ATAC.R)
path_scATAC_remission <- #path of output location of combined processed scATAC data for all patients (Classify_ATAC.R), this could be a normal bone marrow reference as well
path_to_trajectory_peaks <- #path of output location of processed lineage peaks (Analysis_scATAC_Remission_Trajectory_priming.R)
path_to_bgd_peaks <- #path of output location of processed lineage peaks from other lineages (Analysis_scATAC_Remission_Trajectory_priming.R)

#it is recommended to run calculations of priming scores as a new project to prevent overwriting inferred peaks defined during preprocessing

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
require(GenomicRanges)
require(ggnewscale)

# to make the script reproducible
set.seed(42)

scATAC_data<-path_scATAC

trajectory_peaks<-path_to_trajectory_peaks 
bgd_peaks<- path_to_bgd_peaks

scATAC_data<-addPeakSet(ArchRProj = scATAC_data,peakSet = trajectory_peaks,force=T)
scATAC_data<-addPeakMatrix(scATAC_data,force=T)
cell_peaks<-getMatrixFromProject(scATAC_data,"PeakMatrix")
Trajectory_score<-apply(cell_peaks@assays@data$PeakMatrix,2,sum)

scATAC_data<-addPeakSet(ArchRProj = scATAC_data,peakSet = bgd_peaks,force=T)
scATAC_data<-addPeakMatrix(scATAC_data,force=T)
cell_peaks<-getMatrixFromProject(scATAC_data,"PeakMatrix")
bgd_score<-apply(cell_peaks@assays@data$PeakMatrix,2,sum)

lineage_score<-Trajectory_score/bgd_score

scATAC_data<-addCellColData(scATAC_data,data = lineage_score,name = "Trajectory_score",cells =names(scATAC_data),force=T)

#do the same for "normal_cells" for comparison

scATAC_data_HSPC<-path_scATAC_remission

scATAC_data_HSPC<-addCellColData(scATAC_data_HSPC,data = lineage_score,name = "Trajectory_score",cells =names(scATAC_data_HSPC),force=T)

trajectory_peaks<-path_to_trajectory_peaks 
bgd_peaks<- path_to_bgd_peaks

scATAC_data_HSPC<-addPeakSet(ArchRProj = scATAC_data_HSPC,peakSet = trajectory_peaks,force=T)
scATAC_data_HSPC<-addPeakMatrix(scATAC_data_HSPC,force=T)
cell_peaks<-getMatrixFromProject(scATAC_data_HSPC,"PeakMatrix")
Trajectory_score<-apply(cell_peaks@assays@data$PeakMatrix,2,sum)

scATAC_data_HSPC<-addPeakSet(ArchRProj = scATAC_data_HSPC,peakSet = bgd_peaks,force=T)
scATAC_data_HSPC<-addPeakMatrix(scATAC_data_HSPC,force=T)
cell_peaks<-getMatrixFromProject(scATAC_data_HSPC,"PeakMatrix")
bgd_score<-apply(cell_peaks@assays@data$PeakMatrix,2,sum)

lineage_score_HSPC<-Trajectory_score/bgd_score

scATAC_data_HSPC<-addCellColData(scATAC_data_HSPC,data = lineage_score,name = "Trajectory_score",cells =names(scATAC_data_HSPC),force=T)

#plot a priming map as in the paper
#run first for multiple trajectories to generate something like this

Lineage_score_relative<- lineage_score/ median(lineage_score_HSPC)

df_all<-getCellColData(scATAC_data)
df_all$priming_type<-apply(cbind(df_all$Myeloid_score_relative,df_all$Lymphoid_score_relative,df_all$Erythroid_score_relative),1,which.max)
df_all$priming_tscore<-apply(cbind(df_all$Myeloid_score_relative,df_all$Lymphoid_score_relative,df_all$Erythroid_score_relative),1,max)\

pdf("Priming.pdf",width=5,height=5)
ggplot(df_all)+		
    geom_jitter(
        aes(priming_type, y=log2(priming_tscore), colour = log2(priming_tscore)), 
        filter(df_subgroup, priming_type == "1"), 
        size = 1) +
    scale_colour_gradientn(colours = c("#B1B1B2","#1B75BC", "#1B75BC","#1B75BC"),limits=c(0,5),na.value = "grey") +
    labs(colour = "Myeloid") +
    new_scale_colour() +
    geom_jitter(
        aes(priming_type, y=log2(priming_tscore), colour = log2(priming_tscore)), 
        filter(df_subgroup, priming_type == "2"), 
        size = 1) +
    scale_colour_gradientn(colours = c("#B1B2B1","#00A651", "#00A651","#00A651"),limits=c(0,5),na.value = "grey") +
    labs(colour = "Lymphoid") +
    new_scale_colour() +
    geom_jitter(
        aes(priming_type, y=log2(priming_tscore), colour = log2(priming_tscore)), 
        filter(df_subgroup, priming_type == "3"), 
        size = 1) +
    scale_colour_gradientn(colours = c("#B2B1B1","#EE2A7B", "#EE2A7B","#EE2A7B"),limits=c(0,5),na.value = "grey") +
    labs(colour = "Erythroid") +
    coord_polar()+
    facet_wrap(~subgroup)+
    theme_bw()
dev.off()