# NHP_COVID-19_2

<h4>List of genes used to filter the gene expression matrix:</h4>
MT_genes.txt<br>
IG_genes.txt<br>
TR_genes.txt<br>
Y_genes.txt<br>
protein_coding_seurat.txt<br>

<h4>Files for bulk RNA-Seq based AM & non-AM reference used with SingleR</h4>
rlog_IM_AM_AU.txt<br>
pData.txt<br>

<h4>AM and non-AM markers used by SingleR</h4>
AM_nonAM_SingelR_markers.txt<br>
nonAM_AM_SingelR_markers.txt<br>

<h4>Script for bulk RNA-Seq analysis:</h4>
Mmul10-MT246667_Bulk_Analysis.Rmd<br>

<h4>preprocess_BAL.R</h4>
-Reads in filtered feature barcode h5 files for 3 untreated and 2 baricitinib treated BAL samples<br>
-Filters genes - Ig, TCR, mitochondrial, RPS/RPL genes, SARS-CoV2 genes, non-coding<br>
-Creates a list of seurat objects each of which is annotated with SingleR, CellCycleScore, and normalized with SCTransform<br>
-Integrate all the seurat objects using default CCA method<br>
-Subset macrophages and monocytes - split by sample and integrate the macrophages/monocytes again<br>

<h4>Lung_control.R</h4>
-Processes three BAL samples from healthy rhesus (GEO: GSE149758)<br>
-Create seurat objects separately and integrate using CCA<br>
-Manually select macrophage/monocyte cluster<br>
-Subset, split and re-integrate. Discard DC cluster and repeat<br>
-Annotate seurat clusters based on experession of marker genes and use as reference<br>

<h4>BAL_untreated.R</h4>
-Uses the rhesus lung macrophage/monocyte reference to annotate BAL macrophages/monocytes from untreated samples<br>
-Use the bulk RNA-Seq AM and IM reference to classify BAL macrophages/monocytes into AM and non-AM<br>
-Plot UMAP, proportions, percentage reads for pro-inflammatory genes, differential expression<br>

<h4>BAL_treated_untreated.R</h4>
-Uses the rhesus lung macrophage/monocyte reference to annotate BAL macrophages/monocytes from untreated<br>
 and baricitnib treated samples<br>
-Use the bulk RNA-Seq AM and IM reference to classify BAL macrophages/monocytes into AM and non-AM<br>
-Plot UMAP, proportions, etc<br>

<h4>pDC.R</h4>
-Annotate pDC and plot absolute counts and proportions<br>

<h4>mDC.R</h4>
-Annotate mDC, plot absolute counts and proportions, expression of cytokines, chemokines<br>

<h4>Human_Lung_reference.R</h4>
-Use rds object from GEO: GSE135893<br>
-Subset macrophage/monocyte from Control samples and filter some samples<br>
-Split and reintegrate macrophage/monocyte<br>
-Annotate seurat clusters based on expression of marker genes and use as reference<br>

<h4>Healthy_rhesus_human.R</h4>
-Integrate healthy lung macrophage/monocyte datasets from both rhesus and human<br>

<h4>preprocess_Human_BAL.R</h4>
-Preprocess samples from GSE145926<br>
-Integrate all samples and subset the macrophages/monocytes <br>

<h4>Map_HumanBAL_HumanLung.R</h4>
-Use healthy lung macrophage/monocyte reference to annotate human BAL samples<br>
-Plot UMAP, proprotions, percentage of reads, etc.<br>




