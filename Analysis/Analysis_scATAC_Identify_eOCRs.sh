#to generate eOCRs (bash)
#run Classify_ATAC before using

#require before starting (within R): 
#A position sorted bamfile 
#A list of barcodes to include i.e. all barcodes identified as HSPCs in a patient you can subset getCellColdata(scATAC_data) for this
#A list of peaks identified in the data using MACS2 we recommend subsetting the ArchR project to contain a similar number of cells i.e. 1000 of each group you want to compare
#i.e.Use getPeakSet(scATAC_data) for this after identifying the regions:
#scATAC_data <- addReproduciblePeakSet(ArchRProj = scATAC_data, groupBy = "Cells_of_interest", pathToMacs2 = MACS2_location, minCells = 50,method = "p",cutOff = 1e-5,peakMethod = "Macs2")
#then select peaks that are not promoters (getPeakSet(scATAC_data)$peakType != "Promoter") 

#inputfiles
inputbam= cellranger_output/possorted_bam.bam #standard cellranger outputfile containing all alignments
barcode_list= list_of_barcodes_to_call_eOCRs_in.tsv # a list of valid barcodes to get from the cellranger alignment
peak_set= peakset.gff #a list of peaks derived from scATAC data
outputbam= pseudobulk.bam #the outputfile containing pseudobulk bamfiles of selected cells
reads=50000000 # 50 mil reads are recommended for a proper atac analysis, increase this number if the subsetted bams are much larger in the comparison
rose_outdir= output_directory # Rose outputs a lot of files so a new directory is advised

#subset a bamfile to contain only barcodes of interest
subset-bam_linux --bam inputbam --cell-barcodes barcode_list --out-bam pseudobulk.bam --cores 10

#downsample bam files to the same depth
fraction=$(samtools idxstats pseudobulk.bam | cut -f3 | awk -v ct=$reads 'BEGIN {total=0} {total += $1} END {print ct/total}')
samtools view -b -s ${fraction} pseudobulk.bam > pseudobulk.bam_subsampled.bam

#run rose2 
rose2 -i peakset.gff -r pseudobulk.bam_subsampled.bam -g HG38 -o $outdir -t 2500

#within R
#use reduce within GenomicRanges to generate a set of eOCRs to compare across samples