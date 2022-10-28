library("Seurat")
library("ggplot2")
library("dplyr")

setwd("/data/Analysis_BAL_Lung_Human_2022")

options(future.globals.maxSize= 33554432000) # 32000*1024^2

rhesus_ref <- readRDS("/data/Lungs_int3.rds")
colnames(rhesus_ref@meta.data)

DefaultAssay(rhesus_ref) <- "RNA"
list_rhesus <- SplitObject(rhesus_ref, split.by = "Sample")
for (i in names(list_rhesus)) {
  list_rhesus[[i]]  <- DietSeurat(list_rhesus[[i]], counts = TRUE, data = TRUE, scale.data = FALSE, assays = "RNA")
  list_rhesus[[i]] <- SCTransform(list_rhesus[[i]], verbose = FALSE)
}

human_ref <- readRDS("/data/lung_int_kropski_VUH.rds")
colnames(human_ref@meta.data)
DefaultAssay(human_ref) <- "RNA"
list_human <- SplitObject(human_ref, split.by = "Sample_Name")
for (i in names(list_human)) {
  list_human[[i]]  <- DietSeurat(list_human[[i]], counts = TRUE, data = TRUE, scale.data = FALSE, assays = "RNA")
  list_human[[i]] <- SCTransform(list_human[[i]], verbose = FALSE)
}

list_human_rhesus <- c(list_human,list_rhesus)
names(list_human_rhesus)
################################################################################################################################################################

# Reference based integration - human reference

human_ref_int.features <- SelectIntegrationFeatures(object.list = list_human_rhesus, normalization.method = "SCT",nfeatures = 3000, verbose = F) 
human_ref_int.list <- PrepSCTIntegration(list_human_rhesus, anchor.features = human_ref_int.features, verbose = F)
human_ref_int.list <- lapply(X = human_ref_int.list, FUN = RunPCA, features = human_ref_int.features)
anchors <- FindIntegrationAnchors(human_ref_int.list, normalization.method = "SCT",anchor.features = human_ref_int.features,
                                  reduction = "rpca", reference = c(1:6), verbose = F)
human_ref_int <- IntegrateData(anchorset = anchors, normalization.method = "SCT", verbose = F)
DefaultAssay(human_ref_int) <- "integrated"
human_ref_int <- RunPCA(human_ref_int, verbose = FALSE)
human_ref_int <- RunUMAP(human_ref_int, dims = 1:30, verbose = FALSE)
human_ref_int <- FindNeighbors(human_ref_int, dims = 1:30, verbose = FALSE)
human_ref_int <- FindClusters(human_ref_int, verbose = FALSE, resolution = 0.1)
DimPlot(human_ref_int) & NoAxes()

human_ref_int$celltype <- ifelse(human_ref_int$celltype=="Macrophages",human_ref_int$Type,human_ref_int$celltype)
human_ref_int$celltype <- ifelse(human_ref_int$celltype=="Monocytes",human_ref_int$Type,human_ref_int$celltype)
human_ref_int$celltype <- ifelse(human_ref_int$celltype=="Proliferating Macrophages",human_ref_int$Type,human_ref_int$celltype)

human_ref_int$celltype <- factor(human_ref_int$celltype,levels = c("CD163+MRC1+\nMac","CD163+MRC1+\nTREM2+ Mac","CD163+MRC1-\nMac","CD16+\nMono","FABP4hi","SPP1hi","FCN1hi","Proliferating"))
DimPlot(human_ref_int, split.by = "celltype",ncol = 4) & NoAxes()

saveRDS(human_ref_int,"human_ref_int.rds")
colnames(human_ref_int@meta.data)

unique(human_ref_int$orig.ident)

human_ref_rhesussubset <- subset(human_ref_int, cells = row.names(human_ref_int@meta.data[grepl("Lung_NA",human_ref_int$orig.ident),]))
human_ref_humansubset <- subset(human_ref_int, cells = row.names(human_ref_int@meta.data[!grepl("Lung_NA",human_ref_int$orig.ident),]))

dim(human_ref_rhesussubset)
dim(human_ref_humansubset)

prop.table(table(human_ref_humansubset$seurat_clusters,human_ref_humansubset$celltype),1)*100
prop.table(table(human_ref_rhesussubset$seurat_clusters,human_ref_rhesussubset$celltype),1)*100

saveRDS(human_ref_rhesussubset,"human_ref_rhesussubset.rds")
saveRDS(human_ref_humansubset,"human_ref_humansubset.rds")

################################################################################################################################################################

# Reference based integration - human reference (one-to-one orthologs)

orth <- read.table("/data/Ensembl107_Human_Macaque_Orthology.txt", sep = "\t", header = T)
orth <- orth[,c(1,5,6,7,8,9)]
names(orth) <- c("HumanID","RhesusName","RhesusID","Type","Confidence","HumanName")

orth_121 <- unique(orth[orth$Type == "ortholog_one2one" & orth$Confidence == 1,])

same_names <- orth_121[orth_121$RhesusName == orth_121$HumanName,]

human_genes = row.names(human_ref)
rhesus_genes <- row.names(rhesus_ref)
common <- intersect(human_genes,rhesus_genes)

common_same <- common[common %in% same_names$RhesusName]

list_rhesus_orth <- list()

DefaultAssay(rhesus_ref) <- "RNA"
list_rhesus <- SplitObject(rhesus_ref, split.by = "Sample")
for (i in names(list_rhesus)) {
  DefaultAssay(list_rhesus[[i]]) <- "RNA"
  counts <- GetAssayData(list_rhesus[[i]], assay = "RNA", slot = "counts")
  print(dim(counts))
  keep <- row.names(counts) %in% common_same
  new_obj <- CreateSeuratObject(counts = counts[keep,])
  new_obj <- AddMetaData(new_obj, list_rhesus[[i]]@meta.data)
  new_obj <- AddMetaData(new_obj, list_rhesus[[i]]@meta.data)
  print(dim(new_obj))
  new_obj <- SCTransform(new_obj, verbose = FALSE)
  list_rhesus_orth[[i]] <- new_obj
}

list_human_orth <- list()

DefaultAssay(human_ref) <- "RNA"
list_human <- SplitObject(human_ref, split.by = "Sample_Name")
for (i in names(list_human)) {
  counts <- GetAssayData(list_human[[i]], assay = "RNA", slot = "counts")
  keep <- row.names(counts) %in% common_same
  new_obj <- CreateSeuratObject(counts = counts[keep,])
  new_obj <- AddMetaData(new_obj, list_human[[i]]@meta.data)
  new_obj <- AddMetaData(new_obj, list_human[[i]]@meta.data)
  print(dim(new_obj))
  new_obj <- SCTransform(new_obj, verbose = FALSE)
  list_human_orth[[i]] <- new_obj
}

list_human_rhesus_orth <- c(list_human_orth,list_rhesus_orth)

options(future.globals.maxSize= 33554432000) # 32000*1024^2

human_ref_int_orth.features <- SelectIntegrationFeatures(object.list = list_human_rhesus_orth, normalization.method = "SCT",nfeatures = 3000, verbose = F) 
human_ref_int_orth.list <- PrepSCTIntegration(list_human_rhesus_orth, anchor.features = human_ref_int_orth.features, verbose = F)
human_ref_int_orth.list <- lapply(X = human_ref_int_orth.list, FUN = RunPCA, features = human_ref_int_orth.features)
anchors <- FindIntegrationAnchors(human_ref_int_orth.list, normalization.method = "SCT",anchor.features = human_ref_int_orth.features,
                                  reduction = "rpca", reference = c(1:6), verbose = F)
human_ref_int_orth <- IntegrateData(anchorset = anchors, normalization.method = "SCT", verbose = F)
DefaultAssay(human_ref_int_orth) <- "integrated"
human_ref_int_orth <- RunPCA(human_ref_int_orth, verbose = FALSE)
human_ref_int_orth <- RunUMAP(human_ref_int_orth, dims = 1:30, verbose = FALSE)
human_ref_int_orth <- FindNeighbors(human_ref_int_orth, dims = 1:30, verbose = FALSE)
human_ref_int_orth <- FindClusters(human_ref_int_orth, verbose = FALSE, resolution = 0.1)
DimPlot(human_ref_int_orth) & NoAxes()

human_ref_int_orth$celltype <- ifelse(human_ref_int_orth$celltype=="Macrophages",human_ref_int_orth$Type,human_ref_int_orth$celltype)
human_ref_int_orth$celltype <- ifelse(human_ref_int_orth$celltype=="Monocytes",human_ref_int_orth$Type,human_ref_int_orth$celltype)
human_ref_int_orth$celltype <- ifelse(human_ref_int_orth$celltype=="Proliferating Macrophages",human_ref_int_orth$Type,human_ref_int_orth$celltype)

human_ref_int_orth$celltype <- factor(human_ref_int_orth$celltype,levels = c("CD163+MRC1+\nMac","CD163+MRC1+\nTREM2+ Mac","CD163+MRC1-\nMac","CD16+\nMono","FABP4hi","SPP1hi","FCN1hi","Proliferating"))

human_ref_rhesussubset_orth <- subset(human_ref_int_orth, cells = row.names(human_ref_int_orth@meta.data[grepl("Lung_NA",human_ref_int_orth$orig.ident),]))
human_ref_humansubset_orth <- subset(human_ref_int_orth, cells = row.names(human_ref_int_orth@meta.data[!grepl("Lung_NA",human_ref_int_orth$orig.ident),]))

p5 <- DimPlot(human_ref_int_orth, group.by = "seurat_clusters", label = T, raster = T, size = 5) & NoAxes() & NoLegend()
p6 <- DimPlot(human_ref_int_orth, group.by = "seurat_clusters", split.by = "celltype",
              ncol = 4, raster = T, pt.size = 2) & NoAxes()

p5|p6

dhuman_human <- as.data.frame(prop.table(table(human_ref_humansubset_orth$seurat_clusters,human_ref_humansubset_orth$celltype)[,c(5:8)],1)*100)
dhuman_human$Species <- "Human"
dhuman_human

dhuman_rhesus <- as.data.frame(prop.table(table(human_ref_rhesussubset_orth$seurat_clusters,human_ref_rhesussubset_orth$celltype)[,c(1:4)],1)*100)
dhuman_rhesus$Species <- "Rhesus\nmonkey"
dhuman_rhesus

d <- rbind(dhuman_human,dhuman_rhesus)
names(d ) <- c("Cluster","CellType","Percentage","Species")
d$Cluster <- as.character(d$Cluster)
d$Cluster <- factor(d$Cluster, levels = c("0","3","1","2","4"))

d$Cluster_Species <- paste0(d$Species,"\n",d$Cluster)

ggplot(d,aes(x=Species,y=Percentage,group=Species, fill = CellType)) + 
  geom_bar(stat = "identity",position = "stack", width = 0.55) + 
  scale_fill_manual(values = c("#74a9cf","#d7301f","#fcbba1","#238b45","#045a8d","#a50f15","#fb6a4a","#fec44f")) +
  theme_bw() +
  ylab("Percentage of cell\nsubtypes in cluster") +
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank()) + 
  theme(text = element_text(size=14),
        axis.text.x = element_text(size=12, color = "black"),
        axis.text.y = element_text(size=12, color = "black"),
        strip.background =element_rect(fill="white"))  + 
  facet_wrap(~Cluster, ncol = 5)


################################################################################################################################################################

# UCell

DefaultAssay(human_ref) <- "RNA"
human_markers <- FindAllMarkers(human_ref, test.use = "MAST",only.pos = T) 

DefaultAssay(rhesus_ref) <- "RNA"
rhesus_markers <- FindAllMarkers(rhesus_ref, test.use = "MAST",only.pos = T) 

human_markers_sig <- human_markers[human_markers$p_val_adj < 0.05 & human_markers$avg_log2FC >= log2(1.5),]
rhesus_markers_sig <- rhesus_markers[rhesus_markers$p_val_adj < 0.05 & rhesus_markers$avg_log2FC >= log2(1.5),]

fabp4_markers <- row.names(human_markers_sig[human_markers_sig$cluster == "FABP4hi",])
spp1_markers <- row.names(human_markers_sig[human_markers_sig$cluster == "SPP1hi",])
fcn1_markers <- row.names(human_markers_sig[human_markers_sig$cluster == "FCN1hi",])

mrc1neg_markers <- row.names(rhesus_markers_sig[rhesus_markers_sig$cluster=="CD163+MRC1-\nMac",])
mrc1pos_markers <- row.names(rhesus_markers_sig[rhesus_markers_sig$cluster=="CD163+MRC1+\nMac",])
trem2_markers <- row.names(rhesus_markers_sig[rhesus_markers_sig$cluster=="CD163+MRC1+\nTREM2+ Mac",])

fabp4_markers_common <- fabp4_markers[fabp4_markers %in% common_same]
fcn1_markers_common <- fcn1_markers[fcn1_markers %in% common_same]
spp1_markers_common <- spp1_markers[spp1_markers %in% common_same]

mrc1pos_markers_common <- mrc1pos_markers[mrc1pos_markers %in% common_same]
mrc1neg_markers_common <- mrc1neg_markers[mrc1neg_markers %in% common_same]
trem2_markers_common <- trem2_markers[trem2_markers %in% common_same]

list_common_markers <- list(fabp4_markers_common,fcn1_markers_common,spp1_markers_common,mrc1pos_markers_common,mrc1neg_markers_common,trem2_markers_common)

human_list <- list("FABP4hi" = fabp4_markers_common, "SPP1hi" = spp1_markers_common, "FCN1hi" = fcn1_markers_common)
names(human_list)

rhesus_list <- list("MRC1pos" = mrc1pos_markers_common, "TREM2pos" = trem2_markers_common, "MRC1neg" = mrc1neg_markers_common)
names(rhesus_list)

human_ref <- AddModuleScore_UCell(human_ref, features = rhesus_list)
signature.names.rhesus <- paste0(names(rhesus_list), "_UCell")
p1 <- VlnPlot(human_ref, features = signature.names.rhesus, group.by = "Type",cols = c("#74a9cf","#d7301f","#fcbba1","#238b45"))

rhesus_ref <- AddModuleScore_UCell(rhesus_ref, features = human_list)
signature.names.human <- paste0(names(human_list), "_UCell")
rhesus_ref$celltype <- factor(rhesus_ref$celltype, levels = c(levels = c("CD163+MRC1+\nMac","CD163+MRC1+\nTREM2+ Mac","CD163+MRC1-\nMac","CD16+\nMono")))
p2 <- VlnPlot(rhesus_ref, features = signature.names.human, group.by = "celltype",cols = c("#045a8d","#a50f15","#fb6a4a","#fec44f"))

p1/p2

################################################################################################################################################################


library(patchwork)

# Plots

library(ggplot2)

human_ref_int <- readRDS("human_ref_int.rds")
human_ref_rhesussubset <- readRDS("human_ref_rhesussubset.rds")
human_ref_humansubset <- readRDS("human_ref_humansubset.rds")

human_ref$Type <- factor(human_ref$Type, levels = c("FABP4hi","SPP1hi","FCN1hi","Proliferating"))
Idents(human_ref) <- "Type"

DefaultAssay(rhesus_ref) <- "RNA" 
DefaultAssay(human_ref) <- "RNA" 

# Fig 5a & 5b

p1 <- DimPlot(rhesus_ref, group.by =  "celltype", repel = TRUE, 
        cols = c("#045a8d","#a50f15","#fb6a4a","#fec44f"),pt.size = 1, raster=T) + ggtitle("Healthy lung - rhesus macaque")
p2 <- DimPlot(human_ref, group.by =  "Type", repel = TRUE,
              cols = c("#74a9cf","#d7301f","#fcbba1","#238b45"),pt.size = 1, raster=T)  + ggtitle("Healthy lung - human")
p2 + p1 & NoAxes()


# Fig 5c
genes2 <- c("CD163","MRC1","MARCO","FABP4","INHBA","TREM2","SPP1","MERTK","LGMN","SIGLEC10","IL1B","FCN1","VCAN","S100A8","S100A9")
p3 <- DotPlot(rhesus_ref, features = genes2, assay = "RNA") + coord_flip()
p4 <- DotPlot(human_ref, features = genes2, assay = "RNA") + coord_flip()
p4 + p3

# Supplementary Fig 7
p5 <- DimPlot(human_ref_int, group.by = "seurat_clusters", label = T, raster = T, size = 5) & NoAxes() & NoLegend()
p6 <- DimPlot(human_ref_int, group.by = "seurat_clusters", split.by = "celltype",
              ncol = 4, raster = T, pt.size = 2) & NoAxes()

prop.table(table(human_ref_humansubset$seurat_clusters,human_ref_humansubset$celltype),1)*100
dhuman_human <- as.data.frame(prop.table(table(human_ref_humansubset$seurat_clusters,human_ref_humansubset$celltype)[,c(5:8)],1)*100)
dhuman_human$Species <- "Human"
dhuman_human

dhuman_rhesus <- as.data.frame(prop.table(table(human_ref_rhesussubset$seurat_clusters,human_ref_rhesussubset$celltype)[,c(1:4)],1)*100)
dhuman_rhesus$Species <- "Rhesus\nmonkey"
dhuman_rhesus

d <- rbind(dhuman_human,dhuman_rhesus)
names(d ) <- c("Cluster","CellType","Percentage","Species")
d$Cluster <- as.character(d$Cluster)
d$Cluster <- factor(d$Cluster, levels = c("0","3","1","2","4"))

d$Cluster_Species <- paste0(d$Species,"\n",d$Cluster)


ggplot(d,aes(x=Species,y=Percentage,group=Species, fill = CellType)) + 
  geom_bar(stat = "identity",position = "stack", width = 0.55) + 
  scale_fill_manual(values = c("#74a9cf","#d7301f","#fcbba1","#238b45","#045a8d","#a50f15","#fb6a4a","#fec44f")) +
  theme_bw() +
  ylab("Percentage of cell\nsubtypes in cluster") +
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank()) + 
  theme(text = element_text(size=14),
        axis.text.x = element_text(size=12, color = "black"),
        axis.text.y = element_text(size=12, color = "black"),
        strip.background =element_rect(fill="white"))  + 
  facet_wrap(~Cluster, ncol = 5)

