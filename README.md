# Scripts used in the paper "A longitudinal single-cell atlas of treatment response in pediatric AML"

Processed data can be found in GEO: GSE235063, GSE235308

These data are split into several files and do not require reprocessing if processed/filtered data are used (Data_processing folder)

RNA:
raw data (all barcodes) and processed data (barcodes after filtering)

genes.tsv (contains the list of genes)

barcodes.tsv (contains the list of barcodes)

matrix.mtx (contains the expression values)

Also contains metadata per cell (for processed barcodes, see below)

matrix file can be imported using your favorite scRNA analysis package, for Seurat use ReadMtx(mtx=matrix, cells=barcodes, features=genes)

Example of metadata:

Cell_Barcode	GEO_ID	Lambo_et_al_ID	Patient_Sample	Library_ID	Counts	Features	Mitochondria_percent	Classified_Celltype	Seurat_Cluster	Malignant	Patient_ID	Biopsy_Origin	Age_Months	Disease_free_days	Clinical_Blast_Percent	Expected_Driving_Aberration	Subgroup	Color_Subgroup	Known_CNVs	Treatment_Outcome	nCount_RNA	nFeature_RNA
AAACCCAAGCGTCAGA-1	PAULSL_Dx_scRNA	AML1	Remission	AML1_DX	6928	1965	8.28521939953811	Monocytes	9	Normal	PAULSL	Marrow	22	2586	NA	KMT2A-PALM2.AKAP2	KMT2A	#8F3FD5	None	Censored	6928	1965

scRNA analysis performed for the paper is heavily reliant on Seurat (https://satijalab.org/seurat/) 
Standard functions included in this package are not included in this repository and can be found using the well annotated tutorials of the Seurat package

ATAC:
filtered barcodes.tsv (contains barcodes after filtering)

filtered fragments.tsv (contains fragments from barcodes after filtering)

filtered metadata.tsv (contains metadata per cell, see below)

peaks.tsv (contains identified peaks using MACS2 based on filtered data)

peak_matrix.mtx (contains insertions within peaks)

unfilteredfragments.tsv (for purposes requiring all the fragments)

Fragments can be converted to arrow files using ArchR using the createArrowFiles function on the fragments, which is compatible with scripts in this repository

Example of metadata:

Cell_Barcode	GEO_ID	Lambo_et_al_ID	Patient_Sample	Library_ID	Total_fragments	ReadsInTSS	ReadsInPeaks	FractionInPeaks	Classified_Celltype	ArchR_clusters	Malignant	Linked_scRNA_barcode	Linked_scRNA_Sample	Linked_scRNA_confidence_score	Patient_ID	Biopsy_Origin	Age_Months	Disease_free_days	Clinical_Blast_Percent	Expected_Driving_Aberration	Subgroup	Color_Subgroup	Known_CNV	Treatment_Outcome
CCTTAATGTCATTGGT-1	PAULSL_Dx_scATAC	AML1	Diagnosis	AML1_DX	99012	49961	114637	0.578927965416936	B.Cell	C7	Malignant	CATACTTAGTACCCTA-1	AML1_REM	0.736947701403939	PAULSL	Marrow	22	2586	NA	KMT2A-PALM2.AKAP2	KMT2A	#8F3FD5	None	Censored

scATAC analysis performed for the paper is heavily reliant on ArchR (https://www.archrproject.com/) 
Standard functions included in this package are not included in this repository and can be found using the well annotated tutorials on the ArchR website