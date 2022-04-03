library("Seurat")
library("ggplot2")
library("dplyr")

options(future.globals.maxSize= 33554432000) # 32000*1024^2

setwd("/data/Analysis_BAL_Lung_Human_2022")

lung_kropski <- readRDS("/data/Human_Lung_Kropski/GSE135893_ILD_annotated_fullsize.rds")
DimPlot(lung_kropski, label = T) & NoLegend()

colnames(lung_kropski@meta.data)
lung_kropski_macro <- subset(lung_kropski, cells = row.names(lung_kropski@meta.data[(lung_kropski$celltype == "Macrophages"  | lung_kropski$celltype == "Monocytes" | lung_kropski$celltype == "Proliferating Macrophages") & lung_kropski$Diagnosis == "Control" ,]))

DimPlot(lung_kropski_macro,group.by = "seurat_clusters")
DefaultAssay(lung_kropski_macro) <- "RNA"
lung_kropski_macro <- NormalizeData(lung_kropski_macro)
lung_kropski_macro <- FindVariableFeatures(lung_kropski_macro, selection.method = "vst")
lung_kropski_macro <- ScaleData(lung_kropski_macro)
lung_kropski_macro  <- RunPCA(lung_kropski_macro, verbose = FALSE)
lung_kropski_macro  <- RunUMAP(lung_kropski_macro, dims = 1:10, verbose = FALSE)
lung_kropski_macro  <- FindNeighbors(lung_kropski_macro, dims = 1:10, verbose = FALSE)
lung_kropski_macro  <- FindClusters(lung_kropski_macro, verbose = FALSE, resolution = 0.1)
DimPlot(lung_kropski_macro,group.by = "Sample_Source", raster = T, pt.size = 0.5) & NoAxes()
saveRDS(lung_kropski_macro,"lung_kropski_macro.rds")

list_lung_kropski <- SplitObject(lung_kropski_macro, split.by = "Sample_Name")
for (i in names(list_lung_kropski)) {
  list_lung_kropski[[i]]  <- DietSeurat(list_lung_kropski[[i]], counts = TRUE, data = TRUE, scale.data = FALSE, assays = "RNA")
  list_lung_kropski[[i]] <- SCTransform(list_lung_kropski[[i]], verbose = FALSE)
}

names(list_lung_kropski)

# list_lung_kropski <- readRDS("list_lung_kropski.rds")
# Using only Vanderbilt site to avoid potential batch effects and VUHD69 due to low number of macrophages/monocytes
list_lung_kropski <- list_lung_kropski[grepl("VUH",names(list_lung_kropski))]
list_lung_kropski["VUHD69"] <- NULL
names(list_lung_kropski)

lung_int_kropski.features <- SelectIntegrationFeatures(object.list = list_lung_kropski, normalization.method = "SCT",nfeatures = 3000) 
lung_int_kropski.list <- PrepSCTIntegration(list_lung_kropski, anchor.features = lung_int_kropski.features, verbose = FALSE)
anchors <- FindIntegrationAnchors(lung_int_kropski.list, normalization.method = "SCT",anchor.features = lung_int_kropski.features, verbose = FALSE, dims = 1:15)
lung_int_kropski <- IntegrateData(anchorset = anchors, normalization.method = "SCT", verbose = FALSE)
DefaultAssay(lung_int_kropski) <- "integrated"
lung_int_kropski <- RunPCA(lung_int_kropski, verbose = FALSE)
lung_int_kropski <- RunUMAP(lung_int_kropski, dims = 1:15, verbose = FALSE)
lung_int_kropski <- FindNeighbors(lung_int_kropski, dims = 1:15, verbose = FALSE)
lung_int_kropski <- FindClusters(lung_int_kropski, verbose = FALSE, resolution = 0.1)
DimPlot(lung_int_kropski, group.by = "celltype", split.by = "celltype")
DimPlot(lung_int_kropski)

DefaultAssay(lung_int_kropski) <-"RNA"
lung_int_kropski <- NormalizeData(lung_int_kropski)
lung_int_kropski <- FindVariableFeatures(lung_int_kropski, selection.method = "vst")
lung_int_kropski <- ScaleData(lung_int_kropski)

lung_int_kropski <- readRDS("lung_int_kropski_VUH.rds")

DimPlot(lung_int_kropski,group.by = "celltype", label = T)
DimPlot(lung_int_kropski,group.by = "seurat_clusters", label = T, label.size = 6, raster = T) & NoAxes() & NoLegend()
DotPlot(lung_int_kropski,features = c("CD68","CD14","FCGR3A","CD163","MRC1","MARCO","FABP4","SIGLEC1","SPP1","TREM2","LGMN","FCN1","VCAN","S100A8","MKI67"), assay = "RNA") + coord_flip()

lung_int_kropski$Type <- "NA"
lung_int_kropski$Type[lung_int_kropski$seurat_clusters== 0] <- "FABP4hi"
lung_int_kropski$Type[lung_int_kropski$seurat_clusters== 1] <- "SPP1hi"
lung_int_kropski$Type[lung_int_kropski$seurat_clusters== 2] <- "FCN1hi"
lung_int_kropski$Type[lung_int_kropski$seurat_clusters== 3] <- "Proliferating"
lung_int_kropski$Type <- factor(lung_int_kropski$Type, levels = c("FABP4hi","SPP1hi","FCN1hi","Proliferating"))

DimPlot(lung_int_kropski,group.by = "Type", label = T, label.size = 4, raster = T) & NoAxes() & NoLegend()

Idents(lung_int_kropski) <- "Type"

saveRDS(lung_int_kropski,"lung_int_kropski_VUH.rds")
saveRDS(list_lung_kropski,"list_lung_kropski_VUH.rds")
