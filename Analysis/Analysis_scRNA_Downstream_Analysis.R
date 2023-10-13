# script for scRNA-based downstream analysis
# requires process_reference.R and processed scRNA data (Classify_RNA.R) to work

path_scRNA <- #path of output location of processed scRNA data for each patient (Classify_RNA.R)
path_model_probabilities
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
require(ggplot2)
require(presto)
require(ComplexHeatmap)
require(magick)
require(zoo)

reference<-readRDS(path_reference_out)

#run markerlist detection

max_markers<-30
Markerlist<-list()
for (celltypes_of_interest)
{
Markers<- FindMarkers(reference, ident.1 = celltypes_of_interest, group.by = 'BioClassification',test.use = "MAST", only.pos = TRUE, min.pct = 0.25, logfc.threshold = 1)
Markers<-Markers[order(Markers$avg_log2FC,decreasing=T),]
Markers<-Markers[1:max_markers]
Markerlist[[celltypes]]<-rownames(Markers)
}

#draw a heatmap of marker expression

scRNA_data<- path_scRNA 
class_stats<-path_classification_stats
mt_class<-as.data.frame(class_stats[,c(celltypes_of_interest)])
mt_rna<- scRNA_data[["RNA"]]@data[unlist(Markerlist), ] %>% as.matrix()

#for this example 7 celltypes are taken
#reordering is mainly for visualization purposes, it just ensures that high scoring cells are more to the left in the heatmap to show more of a gradient
#can be skipped completely or adjusted for other numbers of celltypes
recalc_mt<-as.matrix(matrix_df)
recalc_mt[,1]<-recalc_mt[,1]*10000
recalc_mt[,2]<-recalc_mt[,2]*1000
recalc_mt[,3]<-recalc_mt[,3]*100
recalc_mt[,4]<-recalc_mt[,4]*10
recalc_mt[,5]<-recalc_mt[,5]*1
recalc_mt[,6]<-recalc_mt[,6]*0.1
recalc_mt[,7]<-recalc_mt[,7]*0.01

mt_class_ordered<-mt_class[order(apply(recalc_mt,1,sum),decreasing = T),]
reordering<-match(rownames(mt_class_ordered),colnames(scRNA_data)) #can be skipped

#heatmap
mt_rna<-scale(mt_rna)
mt_rna<-mt_rna[,rownames(mt_class_ordered)] #can be skipped

#remove extreme values for better contrast
mt_rna[mt_rna > quantile(mt_rna,0.95)]<- quantile(mt_rna,0.95)
mt_rna[mt_rna < quantile(mt_rna,0.05)]<- quantile(mt_rna,0.05)

celltypeAnn<- scRNA_data@meta.data$Cell_type
celltypeAnn<-celltypeAnn[reordering]#can be skipped

cellAnn <- ComplexHeatmap::HeatmapAnnotation(Biopsy = as.factor(scRNA_data@meta.data$Biopsy), col=list(Biopsy=colour))
cellAnn<-cellAnn[reordering]#can be skipped

malAnn <- ComplexHeatmap::HeatmapAnnotation(malignant = as.factor(scRNA_data@meta.data$Malignant_cluster), col=list(malignant=colour))
malAnn<-malAnn[reordering]#can be skipped

row_anno<-gsub('[[:digit:]]+', '', names(Markerlist))

pdf(heatmap.pdf),width=14,height=10)
ht1<- Heatmap(mt_rna, name = "Expression",  
		column_split =  factor(celltypeAnn,levels=c(unique(scRNA_data@meta.data$Cell_type))),
        row_split = factor(row_anno,levels= unique(row_anno)),
        cluster_columns = FALSE,
        show_column_dend = FALSE,
        cluster_column_slices = FALSE,
        column_title_gp = gpar(fontsize = 8),
        column_gap = unit(0.5, "mm"),
        cluster_rows = FALSE,
        show_row_dend = FALSE,
        col = paletteContinuous("blueYellow"),
        row_names_gp = gpar(fontsize = 4),
        column_title = "",
        top_annotation = HeatmapAnnotation(foo = anno_block(gp = gpar(fill = as.vector(colour)))),
	    bottom_annotation = cellAnn,
        show_column_names = FALSE,
        use_raster = TRUE,
	    raster_by_magick = TRUE,
	    raster_magick_filter="Lanczos2",
	    raster_resize_mat = TRUE,
	    raster_quality = 5)
ht2<- Heatmap(t(as.matrix(mt_class_ordered)), name= "Prediction_Scores",
		cluster_rows=F,
	    show_row_dend = FALSE,
	    col=c("white","black"),
	    bottom_annotation = malAnn,
	    show_column_names = FALSE,
	    height = unit(3, "cm"))
ht_list = ht1 %v% ht2
draw(ht_list)
dev.off()

#plot a volcano per subgroup colored by markers

scRNA_data<-scRNA_data[,scRNA_data$Malignant_cluster == "Malignant"]

reps<-10
sample_cells<-1000

diff_express<-list()
for (rep in 1:reps)
{
scRNA_data_subset<-scRNA_data[,c(sample(which(scRNA_data$Biopsy == "Primary"),1000),sample(which(scRNA_data$Biopsy == "Relapse"),1000))]
diff_express[[rep]]<-FindMarkers(object = scRNA_data_subset,group.by = "Biopsy",ident.1 = "Primary",ident.2 = "Relapse",test.use = "MAST",min.pct = 0,logfc.threshold = 0)
diff_express[[rep]]$diff<- diff_express$pct.1-diff_express$pct.2
diff_express[[rep]]$gene<-rownames(diff_express)
}

diff_express$celltype<-"other"
for(ct in names(Markerlist))
{
    diff_express[Markerlist[,ct],]$celltype<-ct
}
diff_express$celltype<-factor(diff_express$celltype,levels=c(unique(names(Markerlist)),"other"))

diff_express$alpha<-1
diff_express[diff_express$celltype == "other",]$alpha<-0.2

ggplot(diff_express,aes(x=avg_log2FC,y=-log10(p_val_adj),alpha=alpha,color=celltype,size=0.5+abs(diff)))+
    geom_point()+
    geom_point(data = diff_express[diff_express$celltype != "other",],aes(x=avg_log2FC,y=-log10(p_val_adj),alpha=alpha,color=celltype,size=0.5+abs(diff)))+
    scale_colour_manual(values=colors_cells)+
    theme_bw()+
    xlim(-5,5)
