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

