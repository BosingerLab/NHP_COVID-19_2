library("Seurat")
library("ggplot2")
library("dplyr")
library("stringr")
library("SingleR")

ref <- BlueprintEncodeData()

setwd("/data/Analysis_BAL_Lung_Human_2022")

options(future.globals.maxSize= 33554432000) # 32000*1024^2


################################################################################################################################################

# Create seurat objects for Control samples

# Original Paper: gene number between 200 and 6000, UMI count above 1000 and mitochondrial gene percentage below 0.1.

createobj_human <- function(obj,file){
  s <- gsub("_filtered.*","",str_split_fixed(obj,"/",2))[,1]
  counts <- Read10X_h5(paste0("/data/GSE145926_human_BAL/",file))
  obj <- CreateSeuratObject(counts=counts,project=s)
  
  print(s)
  
  obj[["percent.mito"]] <- PercentageFeatureSet(obj, pattern = "^MT-")
  obj$Sample <- s
  p <- VlnPlot(obj, features = c("nFeature_RNA", "nCount_RNA","percent.mito"),pt.size=0)
  out <- paste0(s,"_QC.png")
  ggsave(p,file=out)
  
  obj <- subset(obj, subset = nFeature_RNA > 100 & nFeature_RNA < 3500 & percent.mito < 10)
  p <- VlnPlot(obj, features = c("nFeature_RNA", "nCount_RNA","percent.mito"),pt.size=0)
  out <- paste0(s,"_QC_filter.png")
  ggsave(p,file=out)

  obj <- SCTransform(obj, verbose = FALSE)
  obj <- RunPCA(obj, verbose = FALSE)
  obj <- RunUMAP(obj, dims = 1:30, verbose = FALSE)
  obj <- FindNeighbors(obj, dims = 1:30, verbose = FALSE)
  obj <- FindClusters(obj, verbose = FALSE)

  p <- DimPlot(obj, label=TRUE) + NoLegend()
  out <- paste0(s,"_UMAP.png")
  ggsave(p,file=out)

  normdata <- GetAssayData(obj, assay="RNA", slot = "data")
  common <- intersect(rownames(normdata),rownames(ref))
  normdata <- normdata[common,]
  ref <- ref[common,]

  pred.hpca <- SingleR(test = normdata, ref = ref,
                       labels = ref$label.main, assay.type.ref = "logcounts")

  obj[["SingleR.main"]] <- pred.hpca$labels

  Idents(obj) <- "SingleR.main"
  p <- DimPlot(obj, label=TRUE)
  out <- paste0(s,"_SingleR_UMAP.png")
  ggsave(p,file=out)
   
  return(obj)
}

bal_files <- list.files(path = "/data/GSE145926_human_BAL/")
bal_human_list <- mapply(bal_files, FUN=createobj_human,bal_files)
saveRDS(bal_human_list,"BAL_human_list.rds")

################################################################################################################################################


#  Integrate Control, Mild & Severe samples

lapply(bal_human_list,DefaultAssay)

int.features_human <- SelectIntegrationFeatures(object.list = bal_human_list, normalization.method = "SCT", nfeatures = 3000, verbose = F) 
int.list_human <- PrepSCTIntegration(bal_human_list, anchor.features = int.features_human, verbose = F)
int.list_human <- lapply(X = int.list_human, FUN = RunPCA, features = int.features_human)
anchors <- FindIntegrationAnchors(int.list_human, normalization.method = "SCT",anchor.features = int.features_human,
                                  reduction = "rpca", verbose = F)
int_human <- IntegrateData(anchorset = anchors, normalization.method = "SCT", verbose = F)
DefaultAssay(int_human) <- "integrated"
int_human <- RunPCA(int_human, verbose = FALSE)
int_human <- RunUMAP(int_human, dims = 1:30, verbose = F)
int_human <- FindNeighbors(int_human, dims = 1:30, verbose = FALSE)
int_human <- FindClusters(int_human, verbose = FALSE)
Idents(int_human) <- "seurat_clusters"
DimPlot(int_human, group.by = "SingleR.main", label = TRUE)

int_human$Manual = "NA"
plot <- DimPlot(int_human, group.by = "SingleR.main", label = TRUE)
select.cells <- CellSelector(plot = plot)
int_human@meta.data[select.cells,]$Manual <- "Myeloid"

int_human <- readRDS("BAL_human_Myeloid_NA.rds")
p1 <- DimPlot(int_human, group.by = "SingleR.main", label = TRUE, raster = T, repel = T)
p2 <- DimPlot(int_human, group.by = "Manual", label = TRUE, raster = T)
p1 + p2 & NoAxes() & NoLegend()

saveRDS(int_human,"BAL_human_Myeloid_NA.rds")

################################################################################################################################################
# Subset macrophages and monocytes

#int_human <- readRDS("BAL_human_Myeloid_NA.rds")

bal_human_macro_mono <- subset(int_human, cells = row.names(int_human@meta.data[int_human$Manual == "Myeloid" & (int_human$SingleR.main == "Monocytes" | int_human$SingleR.main == "Macrophages"),]))
dim(bal_human_macro_mono)

bal_human_macro_mono$Group <-"Severe"
bal_human_macro_mono$Group[grepl("C51|C52|C100",bal_human_macro_mono$Sample)] <- "Healthy"
bal_human_macro_mono$Group[grepl("C141|C142|C144",bal_human_macro_mono$Sample)] <- "Mild"

bal_human_macro_mono$Sample <- factor(bal_human_macro_mono$Sample, levels = c("GSM4475048_C51","GSM4475049_C52","GSM4475050_C100","GSM4339769_C141","GSM4339770_C142","GSM4339772_C144","GSM4339771_C143","GSM4339773_C145","GSM4339774_C146","GSM4475051_C148","GSM4475052_C149","GSM4475053_C152"))

saveRDS(bal_human_macro_mono, "bal_human_macro_mono.rds")

# bal_human_macro_mono <- readRDS("bal_human_macro_mono.rds")

DefaultAssay(bal_human_macro_mono) <- "RNA"
list_bal_human <- SplitObject(bal_human_macro_mono, split.by = "Sample")
for (i in names(list_bal_human)) {
  DefaultAssay(list_bal_human[[i]]) <- "RNA"
  list_bal_human[[i]]  <- DietSeurat(list_bal_human[[i]], counts = TRUE, data = TRUE, scale.data = FALSE, assays = "RNA")
  list_bal_human[[i]] <- SCTransform(list_bal_human[[i]], verbose = FALSE)
}

bal_human_int.features <- SelectIntegrationFeatures(object.list = list_bal_human, normalization.method = "SCT",nfeatures = 3000, verbose = F) 
bal_human_int.list <- PrepSCTIntegration(list_bal_human, anchor.features = bal_human_int.features, verbose = F)
anchors <- FindIntegrationAnchors(bal_human_int.list, normalization.method = "SCT",anchor.features = bal_human_int.features, verbose = F)
bal_human_int <- IntegrateData(anchorset = anchors, normalization.method = "SCT", verbose = F)

DefaultAssay(bal_human_int) <- "integrated"
bal_human_int <- RunPCA(bal_human_int, verbose = FALSE)
bal_human_int <- RunUMAP(bal_human_int, dims = 1:30, verbose = FALSE)
bal_human_int <- FindNeighbors(bal_human_int, dims = 1:30, verbose = FALSE)
bal_human_int <- FindClusters(bal_human_int, verbose = FALSE, resolution = 0.4)
DimPlot(bal_human_int, label = TRUE) 

bal_human_int$Sample <- factor(bal_human_int$Sample, levels = c("GSM4475048_C51","GSM4475049_C52","GSM4475050_C100","GSM4339769_C141","GSM4339770_C142","GSM4339772_C144","GSM4339771_C143","GSM4339773_C145","GSM4339774_C146","GSM4475051_C148","GSM4475052_C149","GSM4475053_C152"))

bal_human_int$Group <-"Severe"
bal_human_int$Group[grepl("C51|C52|C100",bal_human_int$Sample)] <- "Healthy"
bal_human_int$Group[grepl("C141|C142|C144",bal_human_int$Sample)] <- "Mild"

DefaultAssay(bal_human_int) <- "RNA"
bal_human_int <- NormalizeData(bal_human_int)
bal_human_int <- FindVariableFeatures(bal_human_int, selection.method = "vst")
bal_human_int <- ScaleData(bal_human_int)

DimPlot(bal_human_int, split.by = "Sample", ncol = 3)

saveRDS(bal_human_int,"BAL_humn_int_all_DietSeurat.rds")
saveRDS(list_bal_human,"list_bal_human_macromono_DietSeurat.rds")
