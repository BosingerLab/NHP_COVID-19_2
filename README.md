# NHP_COVID-19_2

List of genes used to filter the gene expression matrix:
MT_genes.txt
IG_genes.txt
TR_genes.txt
Y_genes.txt
protein_coding_seurat.txt

Files for bulk RNA-Seq based AM & non-AM reference used with SingleR
rlog_IM_AM_AU.txt
pData.txt

AM and non-AM markers used by SingleR
AM_nonAM_SingelR_markers.txt
nonAM_AM_SingelR_markers.txt

Script for bulk RNA-Seq analysis:
Mmul10-MT246667_Bulk_Analysis.Rmd

preprocess_BAL.R
-Reads in filtered feature barcode h5 files for 3 untreated and 2 baricitinib treated BAL samples
-Filters genes - Ig, TCR, mitochondrial, RPS/RPL genes, SARS-CoV2 genes, non-coding
-Creates a list of seurat objects each of which is annotated with SingleR, CellCycleScore, and normalized with SCTransform
-Integrate all the seurat objects using default CCA method
-Subset macrophages and monocytes - split by sample and integrate the macrophages/monocytes again

Lung_control.R
-Processes three BAL samples from healthy rhesus (GEO: GSE149758)
-Create seurat objects separately and integrate using CCA
-Manually select macrophage/monocyte cluster
-Subset, split and re-integrate. Discard DC cluster and repeat
-Annotate seurat clusters based on experession of marker genes and use as reference

BAL_untreated.R
-Uses the rhesus lung macrophage/monocyte reference to annotate BAL macrophages/monocytes from untreated samples
-Use the bulk RNA-Seq AM and IM reference to classify BAL macrophages/monocytes into AM and non-AM
-Plot UMAP, proportions, percentage reads for pro-inflammatory genes, differential expression

BAL_treated_untreated.R
-Uses the rhesus lung macrophage/monocyte reference to annotate BAL macrophages/monocytes from untreated
 and baricitnib treated samples
-Use the bulk RNA-Seq AM and IM reference to classify BAL macrophages/monocytes into AM and non-AM
-Plot UMAP, proportions, etc

pDC.R
-Annotate pDC and plot absolute counts and proportions

mDC.R
-Annotate mDC, plot absolute counts and proportions, expression of cytokines, chemokines

Human_Lung_reference.R
-Use rds object from GEO: GSE135893
-Subset macrophage/monocyte from Control samples and filter some samples
-Split and reintegrate macrophage/monocyte
-Annotate seurat clusters based on expression of marker genes and use as reference

Healthy_rhesus_human.R
-Integrate healthy lung macrophage/monocyte datasets from both rhesus and human

preprocess_Human_BAL.R
-Preprocess samples from GSE145926
-Integrate all samples and subset the macrophages/monocytes 

Map_HumanBAL_HumanLung.R
-Use healthy lung macrophage/monocyte reference to annotate human BAL samples
-Plot UMAP, proprotions, percentage of reads, etc.




