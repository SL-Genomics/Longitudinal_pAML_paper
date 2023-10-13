# script to analyze eOCRs within HSPCs
# requires processed scATAC data (Classify_ATAC.R) and processed inferred eOCRs to work 

path_scATAC <- #path of output location of processed scATAC data for a comparison (i.e. HSPCs) (Classify_ATAC.R)
path_eOCRs <- #path of output location of processed eOCRs(Analysis_scATAC_Identify_eOCR.sh)

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
require(BSgenome.Hsapiens.UCSC.hg38)
require(matrixStats)

scATAC_data<-path_scATAC
eOCR_list<-path_eOCRs

#replace peaks with eOCRs, It is recommended to do this in a new project
scATAC_data<-addPeakSet(ArchRProj = scATAC_data,peakSet = eOCR_list,force=T)
scATAC_data<-addPeakMatrix(scATAC_data,force=T)

#link peaks to genes and generate a heatmap
aggregate_cells<-1000 #these calculations take long so cells have to be aggregated

scATAC_data <- addPeak2GeneLinks(ArchRProj = scATAC_data,reducedDims = "IterativeLSI",knnIteration = metacells,k = 100)
ATAC_matrix<-assay(readRDS(metadata(getPeak2GeneLinks(scATAC_data, returnLoops = FALSE))$seATAC))#grab the ATAC data, this contains insertions per eOCR

#annotate the aggregate data for heatmap plotting
annotation_matrix<-matrix(ncol=number_of_annotations,nrow=(metacells-1)) #numbe rof annotations is the amount of annotations you want to add

for (acell in 1:(aggregate_cells-1))
{
cell_barcodes<-unlist(readRDS(metadata(getPeak2GeneLinks(scATAC_data, returnLoops = FALSE))$seATAC)@metadata$KNNList[[acell]])
annotation_matrix[acell,1]<-paste0("K_",acell) #used to sort the matrix
annotation_matrix[acell,2]<-names(sort(table(getCellColData(scATAC_data)[cell_barcodes,]$Biopsy),decreasing=T)[1]) #take the most common annotation i.e. Relapse
annotation_matrix[acell,3]<-names(sort(table(getCellColData(scATAC_data)[cell_barcodes,]$Cell_type),decreasing=T)[1]) #take the most common annotation i.e. HSC
#repeat annotations as desired
}

# select eOCRs linked to a nearby gene with a decent correlation, not used for heatmap plotting but useful for downstream analysis
cor_cutoff= 0.3 #lower range of correlations to include in genes linked to peaks
p2geneDF <- metadata(scATAC_data@peakSet)$Peak2GeneLinks
p2geneDF$genename<-mcols(metadata(p2geneDF)$geneSet)$name[p2geneDF$idxRNA]
#link genes to peaks
p2geneDFcor<-p2geneDF[!is.na(p2geneDF$Correlation) & p2geneDF$Correlation >= cor_cutoff,]

#select the best correlated gene in case more link to one peak
for (enhancer_IDX in 1: length(eOCR_list))
	{
	if (nrow(p2geneDFcor[p2geneDFcor$idxATAC == enhancer_IDX,])	>= 1)
	{
	IDX_enhancer<-p2geneDFcor[p2geneDFcor$idxATAC == enhancer_IDX,]
	p2g_detected_Gene[[enhancer_IDX]]<-IDX_enhancer$genename	
	p2g_top_gene[[enhancer_IDX]]<-IDX_enhancer[which.max(IDX_enhancer$Correlation),]$genename
	}	
	if	(nrow(p2geneDFcor[p2geneDFcor$idxATAC == enhancer_IDX,]) == 0)
	{
	p2g_detected_Gene[[enhancer_IDX]]<-NA
	p2g_top_gene[[enhancer_IDX]]<-NA
	}
}

#use the ArchR plotting function to derive the matrices to plot
plot_eOCR_heatmap_1000_K20_matrices<-plotPeak2GeneHeatmap(ArchRProj = scATAC_data,groupBy = "a variable",corCutOff = cor_cutoff,k = 20,returnMatrices = T)#variable is used for supervised analysis, because the matrices are extracted and reclustered this does not matter

#select the top 1000 most variable correlating eOCRs
p2g_select<-list()
for (p2g_set in unique(plot_eOCR_heatmap_1000_K20_matrices$Peak2GeneLinks)$idxATAC)
	{
	p2g_set_select<-plot_eOCR_heatmap_1000_K20_matrices$Peak2GeneLinks[plot_eOCR_heatmap_1000_K20_matrices$Peak2GeneLinks$idxATAC == p2g_set,]
	p2g_select[[p2g_set]]<-rownames(p2g_set_select[which.max(p2g_set_select$Correlation),])
	}	
eOCR_heatmap_ATAC<-plot_eOCR_heatmap_1000_K20_matrices$ATAC$matrix[unlist(p2g_select),]
eOCR_heatmap_RNA<-plot_eOCR_heatmap_1000_K20_matrices$RNA$matrix[unlist(p2g_select),]
eOCR_heatmap_meta<-plot_eOCR_heatmap_1000_K20_matrices$Peak2GeneLinks[unlist(p2g_select),]
plot_peaks<-rownames(eOCR_heatmap_meta[order(eOCR_heatmap_meta$VarQATAC,decreasing=T),])[1:1000]

#Calculate z_scores based on selected eOCRs

CalcRowZscores <- function(ht = NULL, zlim = zlim){
  z <- sweep(ht - rowMeans(ht), 1, matrixStats::rowSds(ht),`/`)
    z[z > max] <- zlim
    z[z < min] <- -zlim
  return(z)
}

zlim<-2 #min/maximum z-score to plot
z_score_heatmap<-CalcRowZscores(eOCR_heatmap_ATAC[plot_peaks,],zlim)
z_score_heatmap<-z_score_heatmap[,str_sort(colnames(z_score_heatmap),numeric=T)]
z_score_heatmap_RNA<-CalcRowZscores(eOCR_heatmap_RNA[plot_peaks,],zlim)
z_score_heatmap_RNA<-z_score_heatmap_RNA[,str_sort(colnames(z_score_heatmap_RNA),numeric=T)]


#pre-generate the heatmap to save the row order to use in the scATAC and scRNA heatmap, mainly used to sort annotations
ht_ATAC<-Heatmap(z_score_heatmap,
	show_row_dend = F,
	show_column_names = FALSE,
	show_row_names = FALSE,
	use_raster=T,
	raster_by_magick = T,
	raster_magick_filter="Lanczos2",
    raster_quality = 5)
	
ht_ATAC_ord<-draw(ht_ATAC)

#save trees seperately for more resolution
pdf("annotations_heatmaps_eORC.pdf",width=8,height=6)
Heatmap(z_score_heatmap,
	show_row_dend = T,
	show_column_names = FALSE,
	show_row_names = FALSE,
	rect_gp = gpar(type = "none"))
dev.off()

# use the row order in annotations
re_anno_matrix<-annotation_matrix[match(column_order,annotation_matrix[,1]),]
re_anno_heatmap<-eOCR_heatmap_meta[row_order,]

#define annotations here, i.e. Biopsy or celltype
cellAnn_bottom <- ComplexHeatmap::HeatmapAnnotation(
	Biopsy=factor(re_anno_matrix[,2]), 
	col=list(Biopsy=Colour))
cellAnn_top <- ComplexHeatmap::HeatmapAnnotation(
	Celltype=factor(re_anno_matrix[,3]), 
	col=list(Celltype=Colour))

#plot the heatmap of insertions at eOCRs here
row_order<-rownames(z_score_heatmap[row_order(ht_ATAC_ord),])
column_order<-colnames(z_score_heatmap[,column_order(ht_ATAC_ord)])

z_score_heatmap_ord<-z_score_heatmap[,column_order]
z_score_heatmap_ord<-z_score_heatmap_ord[row_order,]

ht_ATAC2<-Heatmap(z_score_heatmap_ord,
	top_annotation = cellAnn_bottom,
	show_row_dend = F,
	cluster_columns = F,
	cluster_rows = F,
	bottom_annotation = cellAnn_bottom,
	show_column_names = FALSE,
	show_row_names = FALSE,
	use_raster=T,
	col = paletteContinuous("coolwarm"),
	raster_by_magick = T,
	raster_magick_filter="Lanczos2",
    raster_quality = 5)

#plot the heatmap of linked gene expression here

z_score_heatmap_RNA_re_ord<-z_score_heatmap_RNA[,column_order]
z_score_heatmap_RNA_re_ord<-z_score_heatmap_RNA_re_ord[row_order,]

#plot genes of interest at the side i.e. MEIS1 and HOXA9
gene_of_interest_list<-c("MEIS1","HOXA9")
ha = rowAnnotation(genes = anno_mark(at = c(which(eOCR_heatmap_meta[row_order,]$gene %in% gene_of_interest_list)),
                                     labels = eOCR_heatmap_meta[row_order,][eOCR_heatmap_meta[row_order,]$gene %in% gene_of_interest_list,]$gene))

ht_RNA<-Heatmap(z_score_heatmap_RNA_re_ord,
                top_annotation = cellAnn_bottom,
                bottom_annotation = cellAnn_bottom,
                cluster_rows = F,
                cluster_columns = F,
                show_row_dend = F,
                show_column_names = FALSE,
                show_row_names = FALSE,
                right_annotation= ha,
                use_raster=T,
                col = paletteContinuous("blueYellow"),
                raster_by_magick = T,
                raster_magick_filter="Lanczos2",
                raster_quality = 5)


pdf("heatmaps_RNA_ATAC_eORC.pdf",width=8,height=6)
ht_ATAC2+ht_RNA
dev.off()
