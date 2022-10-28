setwd("/data/BAL")

library(Seurat)
library(SingleR)
library(stringr)
library(ggplot2)
library(scran)
library(SingleCellExperiment)
library(dplyr)
library(tidyr)
library(reshape2)


bal_int <- readRDS("BAL_int.rds")
list_bal <- readRDS("list_bal_macromono_sct.rds")

#-----------------------------------------------------------------------------------------------------------

# Reference for mapping

lungs_ref <- readRDS("/data/Lung/Lungs_int3.rds")
DefaultAssay(lungs_ref) <- "integrated"
lungs_ref <- RunUMAP(lungs_ref, dims = 1:30, reduction = "pca", return.model = TRUE)

DefaultAssay(lungs_ref) <- "RNA"
DotPlot(lungs_ref,features = c("CD163","CD14","FCGR3","MRC1","CHIT1","MARCO","FABP4","TREM2","APOBEC3A","S100A8","LYVE1"))

lungs_ref$celltype[lungs_ref$celltype == "CD163+MRC1+\nMac"] <- "CD163+MRC1+ Mac"
lungs_ref$celltype[lungs_ref$celltype == "CD163+MRC1+\nTREM2+ Mac"] <- "CD163+MRC1+TREM2+ Mac"
lungs_ref$celltype[lungs_ref$celltype == "CD163+MRC1-\nMac"] <- "CD163+MRC1- Mac"
lungs_ref$celltype[lungs_ref$celltype == "CD16+\nMono"] <- "CD16+ Mono"
unique(lungs_ref$celltype)


#-----------------------------------------------------------------------------------------------------------

# Map lungs reference


list_query <- list_bal
names(list_query)

bal_anchors <- list()
for (i in 1:length(list_query)) {
  bal_anchors[[i]] <- FindTransferAnchors(
    reference = lungs_ref,
    query = list_query[[i]],
    reference.reduction = "pca", 
    dims = 1:30
  )
}

for (i in 1:length(list_query)) {
  list_query[[i]] <- MapQuery(
    anchorset = bal_anchors[[i]], 
    query = list_query[[i]],
    reference = lungs_ref,
    refdata = list(celltype = "celltype"),
    reference.reduction = "pca",
    reduction.model = "umap"
  )
}

bal_lungref <- merge(list_query[[1]], list_query[2:length(list_query)], merge.dr = "ref.umap")

DefaultAssay(bal_lungref) <- "RNA"
bal_lungref <- NormalizeData(bal_lungref)
bal_lungref <- FindVariableFeatures(bal_lungref, selection.method = "vst")
bal_lungref <- ScaleData(bal_lungref)

#-----------------------------------------------------------------------------------------------------------

# SingleR

rData <- read.table("/data/rlog_IM_AM_AU.txt", header=T, row.names = NULL)
rownames(rData) <- make.names(rData[,1], unique = TRUE)
rData <- rData[,-1]
head(rData)

pData <- read.table("/data/pData.txt", header=T)
pData
bulk <- SummarizedExperiment(assays = list(logcounts=rData), colData = pData)

normdata <- GetAssayData(bal_lungref, assay="RNA", slot = "data")
common <- intersect(row.names(normdata),row.names(bulk))
normdata <- normdata[common,]
bulk <- bulk[common,]

bulk.main <- SingleR(test = normdata, ref = bulk, 
                     labels = bulk$Type, assay.type.ref = "logcounts")

bal_lungref[["Bulk.main"]] <- bulk.main$labels
bal_lungref[["Bulk.main.pruned"]] <- bulk.main$pruned.labels
bal_lungref[["Bulk.main.AMscores"]] <- bulk.main$scores[,1]
bal_lungref[["Bulk.main.nonAMscores"]] <- bulk.main$scores[,2]

#-----------------------------------------------------------------------------------------------------------

# Set up factors

bal_lungref$Group[bal_lungref$Group == "Treated"] <- "+Baricitinib"
bal_lungref$Group <- factor(bal_lungref$Group, levels=c("Untreated","Baricitinib"))

bal_lungref$Sample <-  factor(bal_lungref$Sample, levels = c("BAL_-5dpi_Untreated_RHz12","BAL_4dpi_Untreated_RHz12","BAL_-5dpi_Untreated_RQv9","BAL_4dpi_Untreated_RQv9","BAL_-5dpi_Untreated_5-215","BAL_4dpi_Untreated_5-215","BAL_-5dpi_Treated_RLf10","BAL_4dpi_Treated_RLf10","BAL_-5dpi_Treated_RVf12","BAL_4dpi_Treated_RVf12"))
bal_lungref$Sample <- gsub("Treated","+Baricitinib",bal_lungref$Sample)
bal_lungref$Sample <- gsub("BAL_","",bal_lungref$Sample)
bal_lungref$Sample <- gsub("_"," ",bal_lungref$Sample)
bal_lungref$Sample <- factor(bal_lungref$Sample, levels =c("-5dpi Untreated RHz12","4dpi Untreated RHz12",
                                                           "-5dpi Untreated RQv9","4dpi Untreated RQv9",
                                                           "-5dpi Untreated 5-215","4dpi Untreated 5-215",
                                                           "-5dpi +Baricitinib RVf12","4dpi +Baricitinib RVf12",
                                                           "-5dpi +Baricitinib RLf10","4dpi +Baricitinib RLf10"))

bal_lungref$predicted.celltype <- factor(bal_lungref$predicted.celltype, levels = c("CD163+MRC1+ Mac","CD163+MRC1+TREM2+ Mac","CD163+MRC1- Mac","CD16+ Mono"))


bal_lungref$Time_Group <- paste0(bal_lungref$Timepoint," ",bal_lungref$Group)
bal_lungref$Time_Group <- factor(bal_lungref$Time_Group, 
                                 levels = c("-5dpi Untreated","4dpi Untreated","-5dpi +Baricitinib","4dpi +Baricitinib"))


bal_lungref$celltype_Group_Time <- paste0(bal_lungref$predicted.celltype,"_",bal_lungref$Group,"_",bal_lungref$Timepoint)


bal_lungref$Time_celltype <- paste0(bal_lungref$Timepoint,"_",bal_lungref$predicted.celltype)


bal_lungref$celltype_Group_Time <- paste0(bal_lungref$predicted.celltype,"_",bal_lungref$Group,"_",bal_lungref$Timepoint)
bal_lungref$celltype_Group_Time <- factor(bal_lungref$celltype_Group_Time, levels = c("CD163+MRC1- Mac_Untreated_-5dpi","CD163+MRC1- Mac_Untreated_4dpi","CD163+MRC1- Mac_+Baricitinib_-5dpi","CD163+MRC1- Mac_+Baricitinib_4dpi","CD163+MRC1+TREM2+ Mac_Untreated_-5dpi","CD163+MRC1+TREM2+ Mac_Untreated_4dpi","CD163+MRC1+TREM2+ Mac_+Baricitinib_-5dpi","CD163+MRC1+TREM2+ Mac_+Baricitinib_4dpi","CD163+MRC1+ Mac_Untreated_-5dpi","CD163+MRC1+ Mac_Untreated_4dpi","CD163+MRC1+ Mac_+Baricitinib_-5dpi","CD163+MRC1+ Mac_+Baricitinib_4dpi","CD16+ Mono_Untreated_-5dpi","CD16+ Mono_Untreated_4dpi"))

bal_lungref$Bulk.main_Time_Group <- paste(bal_lungref$Bulk.main,bal_lungref$Time_Group,sep="\n")
bal_lungref$Bulk.main_Time_Group <- factor(bal_lungref$Bulk.main_Time_Group, levels = c("non-AM\n-5dpi Untreated","non-AM\n4dpi Untreated","non-AM\n-5dpi +Baricitinib","non-AM\n4dpi +Baricitinib","AM\n-5dpi Untreated","AM\n4dpi Untreated","AM\n-5dpi +Baricitinib","AM\n4dpi +Baricitinib"))


saveRDS(bal_lungref,"bal_lungref.rds")

#-----------------------------------------------------------------------------------------------------------

# UMAP

# Fig 4A
Idents(bal_lungref) <- "predicted.celltype"
DimPlot(bal_lungref, reduction = "ref.umap", group.by =  "Time_Group", repel = TRUE, 
        cols =c("#4d9221","#c51b7d","#ff7f00","#bebada"),pt.size = 2,raster = T)  & NoAxes()

# Fig 4B
DimPlot(bal_lungref, reduction = "ref.umap", group.by =  "predicted.celltype", repel = TRUE, split.by = "Time_Group",ncol = 2,
        cols = c("#045a8d","#a50f15","#fb6a4a","#fec44f"),pt.size = 2,raster = T)  & NoAxes()

# Fig S5A
Idents(bal_lungref) <- "predicted.celltype"
DimPlot(bal_lungref, reduction = "ref.umap", group.by =  "predicted.celltype", repel = TRUE, split.by = "Sample",ncol = 2,
        cols = c("#045a8d","#a50f15","#fb6a4a","#fec44f"),pt.size = 2,raster = T)  & NoAxes()

# Fig S5B
Idents(bal_lungref) <- "Bulk.main"
DimPlot(bal_lungref, reduction = "ref.umap", group.by =  "Bulk.main", repel = TRUE, split.by = "Sample",ncol = 2,
        cols = c("#80b1d3","#fb9a99"),pt.size = 2,raster = T)  & NoAxes()

#-----------------------------------------------------------------------------------------------------------

# VlnPlots

# Fig 4 d,e,f

DefaultAssay(bal_lungref) <- "RNA"
Idents(bal_lungref) <- "celltype_Group_Time"

VlnPlot(bal_lungref,features = c("TNF","IL6","IL10","IFNB1","IL1B"), cols = c("#99d3d4","#e5aeae","#8bd0eb","#dec49f","#88aee1"), stack = T, flip = T) +geom_jitter(size = 0.25) + theme(axis.text.x = element_text(angle = 45, hjust=1)) & NoLegend()
VlnPlot(bal_lungref,features = c("CCL4L1","CXCL10","CXCL3","CXCL8"), cols = c("#cee8c2","#d6bee2","#90c2a9","#b1bbd7"), stack = T, flip = T) +geom_jitter(size = 0.25) + theme(axis.text.x = element_text(angle = 45, hjust=1)) & NoLegend()
VlnPlot(bal_lungref,features = c("IFI27","ISG15","ISG20","MX2"), cols = c("#b5bfa2","#c3eaea","#e8d3ca","#a2bebd"), stack = T, flip = T) +geom_jitter(size = 0.25) + theme(axis.text.x = element_text(angle = 45, hjust=1)) & NoLegend()

#-----------------------------------------------------------------------------------------------------------

# DotPlots AM vs non-AM

# Fig S5e

Idents(bal_lungref) <- "Bulk.main_Time_Group"
DefaultAssay(bal_lungref) <- "RNA"

DotPlot(bal_lungref,features = c("TNF","IL6","IL10","IFNB1","IL1B")) +RotatedAxis()
DotPlot(bal_lungref,features = c("CCL4L1","CXCL10","CXCL3","CXCL8")) + RotatedAxis()
DotPlot(bal_lungref,features = c("IFI27","ISG15","ISG20","MX2")) + RotatedAxis()

#-----------------------------------------------------------------------------------------------------------

# Proportions

# Lung based reference

# Stacked bar chart - Day
d <- as.data.frame(as.matrix(prop.table(table(bal_lungref$Time_Group,bal_lungref$predicted.celltype),1) * 100))
head(d)
names(d) <- c("Sample","Type","Percentage")
d <- d %>% separate(col=Sample, into = c("Timepoint","Group"), sep=" ", remove = F)
d

d$Group <- factor(d$Group, levels=c("Untreated","+Baricitinib"))

ggplot(d, aes(fill=Type, y=Percentage, x=Timepoint)) + 
  geom_bar(position="stack", stat="identity") +
  theme_bw() + ylab("Percentage of Macrophages\n & Monocytes") + 
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank()) + 
  theme(text = element_text(size=14),
        axis.text.x = element_text(size=12, color = "black"),
        axis.text.y = element_text(size=12, color = "black")) +
  facet_grid(cols = vars(Group), scales = "free_x", space = "free_x") + RotatedAxis() +
  theme(strip.background =element_rect(fill="white")) +
  scale_fill_manual(values= c("#045a8d","#a50f15","#fb6a4a","#fec44f"))

write.table(d,"Proportion_BAL_trt_untrt_4C.txt", quote=F,sep="\t")


# Stacked bar chart - AnimalID

d <- as.data.frame(as.matrix(prop.table(table(bal_lungref$Sample,bal_lungref$predicted.celltype),1) * 100))
head(d)
names(d) <- c("Sample","Type","Percentage")
d <- d %>% separate(col=Sample, into = c("Timepoint","Group","AnimalID"), sep=" ", remove = F)

d$AnimalID <- factor(d$AnimalID, levels=c("RHz12","RQv9","5-215","RLf10","RVf12"))

ggplot(d, aes(fill=Type, y=Percentage, x=Timepoint)) + 
  geom_bar(position="stack", stat="identity") +
  theme_bw() + ylab("Proportion of Macrophages \n& Monocytes of all cells") + 
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank()) + 
  theme(text = element_text(size=14),
        axis.text.x = element_text(size=12, color = "black"),
        axis.text.y = element_text(size=12, color = "black")) +
  facet_grid(cols = vars(AnimalID), scales = "free_x", space = "free_x") + RotatedAxis() +
  theme(strip.background =element_rect(fill="white")) +
  scale_fill_manual(values= c("#045a8d","#a50f15","#fb6a4a","#fec44f"))

write.table(d,"Proportion_BAL_trt_untrt_S5C.txt", quote=F,sep="\t")


# Individual

d <- as.data.frame(as.matrix(prop.table(table(bal_lungref$Sample,bal_lungref$predicted.celltype),1) * 100))
names(d) <- c("Sample","Type","Percentage")
d <- d %>% separate(col=Sample, into = c("Timepoint","Group","AnimalID"), sep=" ", remove = F)

d$Type <- as.character(d$Type)
d$Type[d$Type == "CD163+MRC1+ Mac"] <- "CD163+MRC1+\nMac"
d$Type[d$Type == "CD163+MRC1+TREM2+ Mac"] <- "CD163+MRC1+\nTREM2+ Mac"
d$Type[d$Type == "CD163+MRC1- Mac"] <- "CD163+MRC1-\nMac"
d$Type[d$Type == "CD16+ Mono"] <- "CD16+\nMono"

d$Type <- factor(d$Type, levels = c("CD163+MRC1+\nMac","CD163+MRC1+\nTREM2+ Mac","CD163+MRC1-\nMac","CD16+\nMono"))

d$Group <- factor(d$Group, levels = c("Untreated","+Baricitinib"))

ggplot(d, aes(y=Percentage, x=Timepoint)) + 
  geom_line(aes(group=AnimalID), color="grey") +
  geom_point(size = 4, shape = 21, aes(fill= Type)) +
  theme_bw() +
  ylab("Percentage of Macrophages\n& Monocytes") + xlab("") +
  scale_y_continuous(expand=expansion(mult=c(0.1,0.2))) +
  scale_fill_manual(values = c("#045a8d","#a50f15","#fb6a4a","#fec44f")) +
  stat_summary(fun= median, fun.min = median, fun.max = median,
               geom = "crossbar", width = 0.75, color = "black", size = 0.5) +
  facet_grid(Group~Type, scales = "free_y") +
  theme(strip.background =element_rect(fill="white")) +
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank()) +
  theme(text = element_text(size=14),
        axis.text.x = element_text(color = "black", size = 12),
        axis.text.y = element_text(size=12, color = "black"))



# SingleR based reference

# Stacked bar chart - AnimalID

d <- as.data.frame(as.matrix(prop.table(table(bal_lungref$Sample,bal_lungref$Bulk.main),1) * 100))
head(d)
names(d) <- c("Sample","Type","Percentage")
d <- d %>% separate(col=Sample, into = c("Timepoint","Group","AnimalID"), sep=" ", remove = F)
d$AnimalID <- factor(d$AnimalID, levels=c("RHz12","RQv9","5-215","RLf10","RVf12"))


ggplot(d, aes(fill=Type, y=Percentage, x=Timepoint)) + 
  geom_bar(position="stack", stat="identity") +
  theme_bw() + ylab("Proportion of Macrophages \n& Monocytes") + 
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank()) + 
  theme(text = element_text(size=14),
        axis.text.x = element_text(size=12, color = "black"),
        axis.text.y = element_text(size=12, color = "black")) +
  facet_grid(cols = vars(AnimalID), scales = "free_x", space = "free_x") + RotatedAxis() +
  theme(strip.background =element_rect(fill="white")) +
  scale_fill_manual(values= c("#80b1d3","#fb9a99"))

write.table(d,"Proportion_BAL_trt_untrt_S5D.txt", quote=F,sep="\t")





