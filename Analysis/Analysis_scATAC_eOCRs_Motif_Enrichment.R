# script to analyze motif enrichment within eOCRs 
# requires processed scATAC data (Classify_ATAC.R), processed inferred eOCRs to work and likely something to subset eOCRs by i.e. a heatmap (Analysis_scATAC_eOCRs_Downstream_Analysis.R)

path_scATAC <- #path of output location of processed scATAC data for a comparison (i.e. HSPCs) (Classify_ATAC.R)
path_eOCRs <- #path of output location of processed eOCRs(Analysis_scATAC_Identify_eOCR.sh), GRanges format

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

# Make sure this project contains the original peak set and not the eOCRs as peaks!!!
scATAC_data<-path_scATAC
eOCR_list<-path_eOCRs

#variables are to select eOCRs, for instance a cluster identified in a heatmap or eOCRs found in normal vs tumor cells
eOCRs_target<-eOCR_list[variable]
eOCRs_background<-eOCR_list[not_variable]

#grab these from the ArchR github (https://github.com/GreenleafLab/ArchR) functions should included in newer versions of the package
#this script is heavily based on the functions and analyses performed there
source("ValidationUtils.R")
source("AnnotationPeaks.R")

.getAssay <- function(se = NULL, assayName = NULL){
  .assayNames <- function(se){
    names(SummarizedExperiment::assays(se))
  }
  if(is.null(assayName)){
    o <- SummarizedExperiment::assay(se)
  }else if(assayName %in% .assayNames(se)){
    o <- SummarizedExperiment::assays(se)[[assayName]]
  }else{
    stop(sprintf("assayName '%s' is not in assayNames of se : %s", assayName, paste(.assayNames(se),collapse=", ")))
  }
  return(o)
}

# motif enrichment in eOCRs
#get the peaks and motif matches within the dataset
#standard peaks are taken here to enriche for open chromatin since eOCRs are large and not uniformly accessible
matches<- getMatches(sc_ATAC_data)
ranges<- getPeakSet(sc_ATAC_data)[queryHits(findOverlaps(getPeakSet(sc_ATAC_data),eOCRs_target,],type="within")),]
background<-getPeakSet(sc_ATAC_data)[queryHits(findOverlaps(getPeakSet(sc_ATAC_data),eOCRs_background,type="within"))]
idx <- unique(queryHits(findOverlaps(matches, ranges, ignore.strand=TRUE)))
idx_bg <- unique(queryHits(findOverlaps(matches, background, ignore.strand=TRUE)))

matches <- .getAssay(matches,  grep("matches", names(assays(matches)), value = TRUE, ignore.case = TRUE))
  
#Compute Total overlaps
  match_idx <- matches[idx, ,drop=FALSE]
  match_idx_bg <- matches[idx_bg, ,drop=FALSE]
  matchCompareTotal <- Matrix::colSums(matchCompare)
  matchBackgroundTotal <- Matrix::colSums(matchBackground)

#Create a dataframe to compare proportions
  pOut <- data.frame(
    feature = colnames(matches),
    CompareFrequency = matchCompareTotal,
    nCompare = nrow(matchCompare),
    CompareProportion = matchCompareTotal/nrow(matchCompare),
    BackgroundFrequency = matchBackgroundTotal,
    nBackground = nrow(matchBackground),
    BackgroundProporition = matchBackgroundTotal/nrow(matchBackground)
  )
  
#Calculate Enrichment within the dataframe
  pOut$Enrichment <- pOut$CompareProportion / pOut$BackgroundProporition
  
#Get P-Values with Hyper Geometric Test
  pOut$mlog10p <- lapply(seq_len(nrow(pOut)), function(x){
    p <- -phyper(pOut$CompareFrequency[x] - 1, 
     pOut$BackgroundFrequency[x], 
     pOut$nBackground[x] - pOut$BackgroundFrequency[x],
     pOut$nCompare[x], 
     lower.tail = FALSE, log.p = TRUE)
    return(p/log(10))
  }) %>% unlist %>% round(4)

#pOut contains p-value of being more prominent in peaks versus background peaks
  pOut$mlog10Padj <- pmax(pOut$mlog10p - log10(ncol(pOut)), 0)
  pOut <- pOut[order(pOut$mlog10p, decreasing = TRUE), , drop = FALSE]


