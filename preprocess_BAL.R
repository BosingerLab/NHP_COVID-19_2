setwd("/data/BAL")

library(Seurat)
library(SingleR)
library(stringr)
library(ggplot2)
library(scran)
library(SingleCellExperiment)
library(dplyr)

ref <- BlueprintEncodeData()

#--------------------------------------------------------------------------------------------------------------------------------------------------------------------------

# Get list of genes for filtering

y_genes <- as.character(unlist(read.table("Y_genes.txt"),use.names = F))
mt_genes <-  as.character(unlist(read.table("MT_genes.txt"),use.names = F))
ig_genes <- as.character(unlist(read.table("IG_genes.txt"),use.names = F))
tr_genes <- as.character(unlist(read.table("TR_genes.txt"),use.names = F))
prot_coding <- as.character(unlist(read.table("protein_coding_seurat.txt"),use.names = F))


#--------------------------------------------------------------------------------------------------------------------------------------------------------------------------

# Get list of input files

data_path <- "/data/10X/BAL"  

dirs <- list.dirs(data_path, recursive = FALSE)
dirs <- dirs[grepl("Mmul10-MT246667", dirs)]
dirs <- dirs[1:10]
dirs

#--------------------------------------------------------------------------------------------------------------------------------------------------------------------------

# make list of seurat objects for QC

list_obj <- list()

runQC <- function(path){
  
  animalid <- str_extract(path,"5-215|RHz12|RQv9|RLf10|RVf12")
  timepoint <- ifelse(grepl("Day4",path),"4dpi","-5dpi")
  group <- ifelse(grepl("5-215|RHz12|RQv9", path) ,"Untreated",
                  ifelse(grepl("RLf10|RVf12",path) ,"Treated","NA"))
  tissue <- "BAL"
  
  s <- paste(tissue,timepoint,group,animalid, sep="_")
  print(path)
  print(s)
  counts <- Read10X_h5(paste0(path,"/outs/filtered_feature_bc_matrix.h5"))
  
  obj <- CreateSeuratObject(counts=counts,project=s)
  
  rps_genes <- row.names(obj)[grepl("^RP[S,L]",row.names(obj))]
  cov2_genes <- row.names(obj)[grepl("SARS-CoV2",row.names(obj))]
  obj[["percent.cov2"]] <- PercentageFeatureSet(obj, features = cov2_genes)
  obj[["percent.hbb"]] <- PercentageFeatureSet(obj, pattern = "^HBB")
  obj[["percent.mito"]] <- PercentageFeatureSet(obj, features = mt_genes)
  obj[["percent.ig"]] <- PercentageFeatureSet(obj, features = ig_genes)
  obj[["percent.rps.rpl"]] <- PercentageFeatureSet(obj, features = rps_genes)
  obj$log10GenesPerUMI <- log10(obj$nFeature_RNA) / log10(obj$nCount_RNA)
  
  meta <- obj@meta.data[,c("percent.hbb","percent.mito","percent.rps.rpl","percent.ig","percent.cov2","log10GenesPerUMI")]
  
  is_coding <- row.names(obj) %in% prot_coding
  non_coding <- row.names(obj)[!is_coding]
  
  filt_genes <- c("HBB",y_genes, mt_genes, rps_genes, ig_genes, tr_genes, cov2_genes, non_coding)
  is_filt <- row.names(obj) %in% filt_genes
  
  obj <- CreateSeuratObject(counts=counts[!is_filt,])
  
  obj$Sample <- s
  obj$Group <- group
  obj$Tissue <- tissue
  obj$Timepoint <- timepoint
  
  obj@meta.data <- transform(merge(obj@meta.data, meta, by=0, ln.x=TRUE), row.names=Row.names, Row.names=NULL)
  
  p <- VlnPlot(obj, features = c("nFeature_RNA", "nCount_RNA","percent.mito","percent.hbb","percent.rps.rpl","percent.ig","percent.cov2"))
  out <- paste0(s,"_QC.png")
  ggsave(p,file=out)
  
  list_obj[[s]] <<- obj
}

lapply(dirs,runQC)



#--------------------------------------------------------------------------------------------------------------------------------------------------------------------------

# Run individual samples through seurat pipeline

lapply(list_obj,dim)

createobj <- function(obj){
  obj <- subset(obj, subset= (nFeature_RNA >= 200) & (nFeature_RNA <=4000) &
                  (log10GenesPerUMI >= 0.8) &
                  (percent.hbb < 10) & 
                  (percent.mito < 20) & 
                  (percent.rps.rpl < 30))
  s <- unique(obj$Sample)
  
  p <- VlnPlot(obj, features = c("nFeature_RNA", "nCount_RNA","percent.mito","percent.hbb","percent.rps.rpl","percent.ig"))
  out <- paste0(s,"_QC_filter.png")
  ggsave(p,file=out)
  
  obj <- NormalizeData(obj)
  obj <- FindVariableFeatures(obj, selection.method = "vst")
  obj <- ScaleData(obj)
  
  DefaultAssay(obj) <- "RNA"
  normdata <- GetAssayData(obj, assay="RNA", slot = "data")
  common <- intersect(row.names(normdata),row.names(ref))
  normdata <- normdata[common,]
  ref <- ref[common,]
  
  bp.main <- SingleR(test = normdata, ref = ref,
                     labels = ref$label.main, assay.type.ref = "logcounts")
  
  obj[["BP.main"]] <- bp.main$labels
  obj[["BP.main.pruned"]] <- bp.main$pruned.labels
  
  bp.fine <- SingleR(test = normdata, ref = ref,
                     labels = ref$label.fine, assay.type.ref = "logcounts")
  obj[["BP.fine"]] <- bp.fine$labels
  obj[["BP.fine.pruned"]] <- bp.fine$pruned.labels
  
  obj <- CellCycleScoring(obj, s.features = cc.genes$s.genes, g2m.features = cc.genes$g2m.genes)
  
  obj <- SCTransform(obj, verbose = FALSE)
  obj <- RunPCA(obj, verbose = FALSE)
  obj <- RunUMAP(obj, dims = 1:30, verbose = FALSE)
  obj <- FindNeighbors(obj, dims = 1:30, verbose = FALSE)
  obj <- FindClusters(obj, verbose = FALSE)
  
  p <- DimPlot(obj, label=TRUE) + NoLegend()
  out <- paste0(s,"_UMAP.png")
  ggsave(p,file=out)
  
  Idents(obj) <- "SingleR"
  p <- DimPlot(obj, label=TRUE)
  out <- paste0(s,"_SingleR_UMAP.png")
  ggsave(p,file=out)
  
  return(obj)  
}

list_obj2 <- lapply(list_obj,createobj)
saveRDS(list_obj2,"BAL_list_objects_all.rds")

#--------------------------------------------------------------------------------------------------------------------------------------------------------------------------

# Integrate all samples

options(future.globals.maxSize= 33554432000) # 32000*1024^2

int.features <- SelectIntegrationFeatures(object.list = list_obj2, normalization.method = "SCT",nfeatures = 3000) 
int.list <- PrepSCTIntegration(list_obj2, anchor.features = int.features, verbose = FALSE)
anchors <- FindIntegrationAnchors(int.list, normalization.method = "SCT",anchor.features = int.features, verbose = FALSE)
int <- IntegrateData(anchorset = anchors, normalization.method = "SCT", verbose = FALSE)
int <- RunPCA(int, verbose = FALSE)
int <- RunUMAP(int, dims = 1:30)
int <- FindNeighbors(int, dims = 1:30, verbose = FALSE)
int <- FindClusters(int, verbose = FALSE)
Idents(int) <- "BP.main"

# Manually select macrophages/monocytes

int$Manual = "NA"
plot <- DimPlot(int, group.by = "BP.main", label = TRUE)
select.cells <- CellSelector(plot = plot)
int@meta.data[select.cells,]$Manual <- "Myeloid"
p1 <- DimPlot(int, group.by = "BP.main", label = TRUE)
p2 <- DimPlot(int, group.by = "Manual", label = TRUE)
p1 + p2

DimPlot(int, group.by = "Manual", split.by = "BP.main", ncol = 4)
saveRDS(int,"BAL_int_filt_Myeloid.rds")

int <- readRDS("BAL_int_filt_Myeloid.rds")

VlnPlot(int,features = "percent.cov2")

as.data.frame(table(int[,int$percent.cov2 > 0]$BP.main))

#--------------------------------------------------------------------------------------------------------------------------------------------------------------------------

# Subset macro_mono

bal_macro_mono <- subset(int, cells = row.names(int@meta.data[int$Manual == "Myeloid" & (int$BP.main == "Monocytes" | int$BP.main == "Macrophages"),]))
DefaultAssay(bal_macro_mono) <- "RNA"
dim(bal_macro_mono)
saveRDS(bal_macro_mono, "BAL_macro_mono.rds")

#--------------------------------------------------------------------------------------------------------------------------------------------------------------------------

# Split and integrate all samples again

# bal_macro_mono <- readRDS("BAL_macro_mono.rds")
list_bal <- SplitObject(bal_macro_mono, split.by = "Sample")
for (i in names(list_bal)) {
  DefaultAssay(list_bal[[i]]) <- "RNA"
  list_bal[[i]]  <- DietSeurat(list_bal[[i]], counts = TRUE, data = TRUE, scale.data = FALSE, assays = "RNA")
  list_bal[[i]] <- SCTransform(list_bal[[i]], method = "glmGamPoi", verbose = FALSE)
}
bal_int.features <- SelectIntegrationFeatures(object.list = list_bal, normalization.method = "SCT",nfeatures = 3000) 
bal_int.list <- PrepSCTIntegration(list_bal, anchor.features = bal_int.features, verbose = FALSE)
anchors <- FindIntegrationAnchors(bal_int.list, normalization.method = "SCT",anchor.features = bal_int.features, verbose = FALSE)
bal_int <- IntegrateData(anchorset = anchors, normalization.method = "SCT", verbose = FALSE)
DefaultAssay(bal_int) <- "integrated"
bal_int <- RunPCA(bal_int, verbose = FALSE)
bal_int <- RunUMAP(bal_int, dims = 1:30, verbose = FALSE)
bal_int <- FindNeighbors(bal_int, dims = 1:30, verbose = FALSE)
bal_int <- FindClusters(bal_int, verbose = FALSE, resolution = 0.1)

DimPlot(bal_int,group.by = "BP.main")
DimPlot(bal_int,group.by = "seurat_clusters")

saveRDS(bal_int,"BAL_int.rds")
saveRDS(list_bal,"list_bal_macromono_sct.rds")


bal_int <- readRDS("BAL_int.rds")
list_bal <- readRDS("list_bal_macromono_sct.rds")

#--------------------------------------------------------------------------------------------------------------------------------------------------------------------------


