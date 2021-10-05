setwd("/data/Lung")
options(future.globals.maxSize= 33554432000) # 32000*1024^2


library(Seurat)
library(SingleR)
library(stringr)
library(ggplot2)
library(scran)
library(SingleCellExperiment)
library(dplyr)

ref <- BlueprintEncodeData()

#------------------------------------------------------------------------------------------------------------------------------------

# Adding some mt genes from Mmul8
mt_genes <-  as.character(unlist(read.table("MT_genes.txt"),use.names = F))
mt_genes<- c(mt_genes,"ND6","CO1","ND5","ND2","ND1","ND3","ATP6","COX2","ND4","ATP8","CO3","CYTB","ND4L")

#------------------------------------------------------------------------------------------------------------------------------------

data_path <- "/data/GSE149758_Lung_control/"  

# Get list of files
lung_files <- list.files(data_path, recursive = FALSE)
lung_files <- lung_files[grepl("_control_", lung_files)]

#------------------------------------------------------------------------------------------------------------------------------------

# Create list to save seurat objects
list_obj_lung <- list()

runQC <- function(path){
  
  animalid <- str_extract(path,"31810|LE99|LD09")
  timepoint <- "NA"
  group <- "Control"
  tissue <- "Lung"
  
  s <- paste(tissue,timepoint,group,animalid, sep="_")
  
  counts <- Read10X_h5(paste0("/data/GSE149758_Lung_control/",path))
  
  obj <- CreateSeuratObject(counts=counts,project=s)
  rps_genes <- row.names(obj)[grepl("^RP[S,L]",row.names(obj))]
  obj[["percent.mito"]] <- PercentageFeatureSet(obj, features = row.names(obj)[(row.names(obj) %in% mt_genes)])
  obj[["percent.rps.rpl"]] <- PercentageFeatureSet(obj, features = rps_genes)
  obj$log10GenesPerUMI <- log10(obj$nFeature_RNA) / log10(obj$nCount_RNA)
  
  obj$Sample <- s
  obj$Group <- group
  obj$Tissue <- tissue
  obj$Timepoint <- timepoint
  
  p <- VlnPlot(obj, features = c("nFeature_RNA", "nCount_RNA","percent.mito","percent.rps.rpl"))
  out <- paste0(s,"_QC.png")
  ggsave(p,file=out)
  
  list_obj_lung[[s]] <<- obj
}

# Run QC and create seurat objects
lapply(lung_files,runQC)

#------------------------------------------------------------------------------------------------------------------------------------

# Function to filter cells, run SingleR, run seurat pipeline on individual samples
createobj <- function(obj){
  obj <- subset(obj, subset= (nFeature_RNA >= 200) & (nFeature_RNA <=4000) &
                  (log10GenesPerUMI >= 0.8) &
                  (percent.mito < 20) & 
                  (percent.rps.rpl < 50))
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

list_obj2_lung <- lapply(list_obj_lung,createobj)
saveRDS(list_obj2_lung,"Lung_list_objects_all.rds")

#------------------------------------------------------------------------------------------------------------------------------------

# Integrate Samples

int_lungs.features <- SelectIntegrationFeatures(object.list = list_obj2_lung, normalization.method = "SCT",nfeatures = 3000) 
int_lungs.list <- PrepSCTIntegration(list_obj2_lung, anchor.features = int_lungs.features, verbose = FALSE)
anchors <- FindIntegrationAnchors(int_lungs.list, normalization.method = "SCT",anchor.features = int_lungs.features, verbose = FALSE)
int_lungs <- IntegrateData(anchorset = anchors, normalization.method = "SCT", verbose = FALSE)
int_lungs <- RunPCA(int_lungs, verbose = FALSE)
int_lungs <- RunUMAP(int_lungs, dims = 1:30)
int_lungs <- FindNeighbors(int_lungs, dims = 1:30, verbose = FALSE)
int_lungs <- FindClusters(int_lungs, verbose = FALSE)
Idents(int_lungs) <- "BP.main"
DimPlot(int_lungs, group.by = "BP.main", label = TRUE)

# Manually select macrophage/monocyte cluster

int_lungs$Manual = "NA"
plot <- DimPlot(int_lungs, group.by = "BP.main", label = TRUE)
select.cells <- CellSelector(plot = plot)
int_lungs@meta.data[select.cells,]$Manual <- "Myeloid"
# int_lungs@meta.data[select.cells,]$Manual <- "NA"

DimPlot(int_lungs, group.by = "BP.main", label = TRUE, split.by = "Manual")
DimPlot(int_lungs, group.by = "BP.main", split.by = "Manual" ,label = TRUE)
saveRDS(int_lungs,"Lungs_Myeloid_NA.rds")

#------------------------------------------------------------------------------------------------------------------------------------

# Subset Get macrophages & monocytes
lungs_macro_mono <- subset(int_lungs, cells = row.names(int_lungs@meta.data[int_lungs$Manual == "Myeloid" & (int_lungs$BP.main == "Monocytes" | int_lungs$BP.main == "Macrophages"),]))
dim(lungs_macro_mono)
saveRDS(lungs_macro_mono, "lungs_macro_mono.rds")

#------------------------------------------------------------------------------------------------------------------------------------

# Split and integrate all samples

list_lungs <- SplitObject(lungs_macro_mono, split.by = "Sample")
for (i in names(list_lungs)) {
  DefaultAssay(list_lungs[[i]]) <- "RNA"
  list_lungs[[i]]  <- DietSeurat(list_lungs[[i]], counts = TRUE, data = TRUE, scale.data = FALSE, assays = "RNA")
  list_lungs[[i]] <- SCTransform(list_lungs[[i]], method = "glmGamPoi", verbose = FALSE)
}
lungs_int.features <- SelectIntegrationFeatures(object.list = list_lungs, normalization.method = "SCT",nfeatures = 3000) 
lungs_int.list <- PrepSCTIntegration(list_lungs, anchor.features = lungs_int.features)
lungs_int.anchors <- FindIntegrationAnchors(lungs_int.list, normalization.method = "SCT",anchor.features = lungs_int.features)
lungs_int <- IntegrateData(anchorset = lungs_int.anchors, normalization.method = "SCT")
DefaultAssay(lungs_int) <- "integrated"
lungs_int <- RunPCA(lungs_int, verbose = FALSE)
lungs_int <- RunUMAP(lungs_int, dims = 1:40, verbose = FALSE)
lungs_int <- FindNeighbors(lungs_int, dims = 1:40, verbose = FALSE)
lungs_int <- FindClusters(lungs_int, verbose = FALSE, resolution = 0.1)
DimPlot(lungs_int, label = T)

DefaultAssay(lungs_int) <- "RNA"
lungs_int <- NormalizeData(lungs_int)
lungs_int <- FindVariableFeatures(lungs_int, selection.method = "vst")
lungs_int <- ScaleData(lungs_int)

lung_markers <- FindAllMarkers(lungs_int,test.use = "MAST", only.pos = T, logfc.threshold = 0.25, min.pct = 0.25)

lung_markers %>%
  group_by(cluster) %>%
  top_n(n = 10, wt = avg_log2FC) -> top10
DoHeatmap(lungs_int, features = top10$gene) + NoLegend()

VlnPlot(lungs_int,c("CD163","MRC1","FABP4","TREM2","S100A8","CD14","FCGR3","FCER1A","CD1C","GZMB","FLT3","CLEC9A"))

write.table(lung_markers,"Lung_macromono_round1_markers.txt", quote=F, sep="\t")

saveRDS(list_lungs,"list_lungs_macromono_round1.rds")
saveRDS(lungs_int,"Lungs_int_round1.rds")

#------------------------------------------------------------------------------------------------------------------------------------

# Remove DC cluster and integrate all the sample again


lungs_int2 <- subset(lungs_int,cells = row.names( lungs_int@meta.data[(lungs_int$seurat_clusters == "0" | lungs_int$seurat_clusters == "1" | lungs_int$seurat_clusters == "2" | lungs_int$seurat_clusters == "3"),]))
list_lungs2 <- SplitObject(lungs_int2, split.by = "Sample")

# list_lungs2 <- readRDS("list_lungs_macromono2.rds")

for (i in names(list_lungs2)) {
  DefaultAssay(list_lungs2[[i]]) <- "RNA"
  list_lungs2[[i]]  <- DietSeurat(list_lungs2[[i]], counts = TRUE, data = TRUE, scale.data = FALSE, assays = "RNA")
  list_lungs2[[i]] <- SCTransform(list_lungs2[[i]], method = "glmGamPoi", verbose = FALSE)
}
lungs_int2.features <- SelectIntegrationFeatures(object.list = list_lungs2, normalization.method = "SCT", nfeatures = 3000, verbose = FALSE) 
lungs_int2.list <- PrepSCTIntegration(list_lungs2, anchor.features = lungs_int2.features, verbose = FALSE)
anchors <- FindIntegrationAnchors(lungs_int2.list, normalization.method = "SCT",anchor.features = lungs_int2.features, verbose = FALSE)
lungs_int2 <- IntegrateData(anchorset = anchors, normalization.method = "SCT", verbose = FALSE)
DefaultAssay(lungs_int2) <- "integrated"
lungs_int2 <- RunPCA(lungs_int2, verbose = FALSE)
lungs_int2 <- RunUMAP(lungs_int2, dims = 1:30, verbose = FALSE)
lungs_int2 <- FindNeighbors(lungs_int2, dims = 1:30, verbose = FALSE)
lungs_int2 <- FindClusters(lungs_int2, verbose = FALSE, resolution = 0.1)
Idents(lungs_int2) <- "seurat_clusters"
DimPlot(lungs_int2)

DefaultAssay(lungs_int2) <- "RNA"
lungs_int2 <- NormalizeData(lungs_int2)
lungs_int2 <- FindVariableFeatures(lungs_int2, selection.method = "vst")
lungs_int2 <- ScaleData(lungs_int2)

FeaturePlot(lungs_int2,c("CD163","MRC1","MARCO","TREM2","APOBEC3A","S100A8","FCGR3"))
DotPlot(lungs_int2,features = c("CD163","MRC1","MARCO","TREM2","APOBEC3A","S100A8","FCGR3","CD14","FCER1A","CD1C")) + coord_flip()

lung_markers2 <- FindAllMarkers(lungs_int2,test.use = "MAST", only.pos = T)

lung_markers2 %>%
  group_by(cluster) %>%
  top_n(n = 10, wt = avg_log2FC) -> top10
DoHeatmap(lungs_int2, features = top10$gene) + NoLegend()

saveRDS(list_lungs2,"list_lungs_macromono2.rds")
saveRDS(lungs_int2,"Lungs_int2.rds")


lungs_int3 <- subset(lungs_int2,cells = row.names( lungs_int2@meta.data[(lungs_int2$seurat_clusters != "4" ),]))
list_lungs3 <- SplitObject(lungs_int3, split.by = "Sample")

# list_lungs3 <- readRDS("list_lungs_macromono2.rds")

for (i in names(list_lungs3)) {
  DefaultAssay(list_lungs3[[i]]) <- "RNA"
  list_lungs3[[i]]  <- DietSeurat(list_lungs3[[i]], counts = TRUE, data = TRUE, scale.data = FALSE, assays = "RNA")
  list_lungs3[[i]] <- SCTransform(list_lungs3[[i]], method = "glmGamPoi", verbose = FALSE)
}
lungs_int3.features <- SelectIntegrationFeatures(object.list = list_lungs3, normalization.method = "SCT", nfeatures = 3000, verbose = FALSE) 
lungs_int3.list <- PrepSCTIntegration(list_lungs3, anchor.features = lungs_int3.features, verbose = FALSE)
anchors <- FindIntegrationAnchors(lungs_int3.list, normalization.method = "SCT",anchor.features = lungs_int3.features, verbose = FALSE)
lungs_int3 <- IntegrateData(anchorset = anchors, normalization.method = "SCT", verbose = FALSE)
DefaultAssay(lungs_int3) <- "integrated"
lungs_int3 <- RunPCA(lungs_int3, verbose = FALSE)
lungs_int3 <- RunUMAP(lungs_int3, dims = 1:30, verbose = FALSE)
lungs_int3 <- FindNeighbors(lungs_int3, dims = 1:30, verbose = FALSE)
lungs_int3 <- FindClusters(lungs_int3, verbose = FALSE, resolution = 0.1)
Idents(lungs_int3) <- "seurat_clusters"
DimPlot(lungs_int3)

DefaultAssay(lungs_int3) <- "RNA"
lungs_int3 <- NormalizeData(lungs_int3)
lungs_int3 <- FindVariableFeatures(lungs_int3, selection.method = "vst")
lungs_int3 <- ScaleData(lungs_int3)

FeaturePlot(lungs_int3,c("CD163","MRC1","MARCO","TREM2","APOBEC3A","S100A8","CD14","FCGR3"))
DotPlot(lungs_int3,features = c("CD163","MRC1","MARCO","TREM2","APOBEC3A","S100A8","FCGR3","CD14","SIGLEC1")) + coord_flip()
VlnPlot(lungs_int3,features = c("CD163","MRC1","MARCO","TREM2","APOBEC3A","S100A8","CD14","FCGR3"))

lungs_int3$celltype = "NA"
lungs_int3$celltype[lungs_int3$seurat_clusters == 0] <- "CD163+MRC1+ Mac"
lungs_int3$celltype[lungs_int3$seurat_clusters == 1] <- "CD163+MRC1- Mac"
lungs_int3$celltype[lungs_int3$seurat_clusters == 2] <- "CD163+MRC1+TREM2+ Mac"
lungs_int3$celltype[lungs_int3$seurat_clusters == 3] <- "CD16+ Mono"

Idents(lungs_int3) <- "seurat_clusters"
DimPlot(lungs_int3)

DefaultAssay(lungs_int3) <- "RNA"
FeaturePlot(lungs_int3,c("CD163","CD14","FCGR3","MRC1","CHIT1","MARCO","TREM2","APOBEC3A","S100A8"))
VlnPlot(lungs_int3,features = c("CD163","CD14","FCGR3","MRC1","CHIT1","MARCO","TREM2","APOBEC3A","S100A8"))

Idents(lungs_int3) <- "celltype"
DotPlot(lungs_int3,features = c("CD163","CD14","FCGR3","MRC1","CHIT1","MARCO","TREM2","APOBEC3A","S100A8")) 


# SingleR

rData <- read.table("/data/rlog_IM_AM_AU.txt", header=T, row.names = NULL)
rownames(rData) <- make.names(rData[,1], unique = TRUE)
rData <- rData[,-1]
head(rData)

pData <- read.table("/data/pData.txt", header=T)
pData
bulk <- SummarizedExperiment(assays = list(logcounts=rData), colData = pData)

normdata <- GetAssayData(lungs_int3, assay="RNA", slot = "data")
common <- intersect(row.names(normdata),row.names(bulk))
normdata <- normdata[common,]
bulk <- bulk[common,]

bulk.main <- SingleR(test = normdata, ref = bulk, 
                     labels = bulk$Type, assay.type.ref = "logcounts")

lungs_int3[["Bulk.main"]] <- bulk.main$labels
lungs_int3[["Bulk.main.pruned"]] <- bulk.main$pruned.labels

p1 <- DimPlot(lungs_int3, group.by = "celltype", ncol = 1)
p2 <- DimPlot(lungs_int3, group.by = "Bulk.main", ncol = 1)
p3 <- DimPlot(lungs_int3, group.by = "Bulk.main", split.by = "celltype", ncol = 2)

saveRDS(list_lungs3,"list_lungs_macromono3.rds")
saveRDS(lungs_int3,"Lungs_int3.rds")

lungs_int3 <- readRDS("Lungs_int3.rds")
