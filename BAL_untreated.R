setwd("/yerkes-cifs/runs/Analysis/COVID_NHP_Study1/R4_10X_Seuratv4_v2/BAL")

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

lungs_ref <- readRDS("/yerkes-cifs/runs/Analysis/COVID_NHP_Study1/R4_10X_Seuratv4_v2/Lung/Lungs_int3.rds")
DefaultAssay(lungs_ref) <- "integrated"
lungs_ref <- RunUMAP(lungs_ref, dims = 1:30, reduction = "pca", return.model = TRUE)

lungs_ref$celltype[lungs_ref$celltype == "CD163+MRC1+\nMac"] <- "CD163+MRC1+ Mac"
lungs_ref$celltype[lungs_ref$celltype == "CD163+MRC1+\nTREM2+ Mac"] <- "CD163+MRC1+TREM2+ Mac"
lungs_ref$celltype[lungs_ref$celltype == "CD163+MRC1-\nMac"] <- "CD163+MRC1- Mac"
lungs_ref$celltype[lungs_ref$celltype == "CD16+\nMono"] <- "CD16+ Mono"

DefaultAssay(lungs_ref) <- "integrated"
lungs_ref <- RunUMAP(lungs_ref, dims = 1:30, reduction = "pca", return.model = TRUE)


#-----------------------------------------------------------------------------------------------------------

# Map lungs reference

list_unt_query <- list_query[grepl("Untreated",names(list_bal))]
names(list_unt_query)

bal_unt_anchors <- list()
for (i in 1:length(list_unt_query)) {
  bal_unt_anchors[[i]] <- FindTransferAnchors(
    reference = lungs_ref,
    query = list_unt_query[[i]],
    reference.reduction = "pca", 
    dims = 1:30
  )
}

for (i in 1:length(list_unt_query)) {
  list_unt_query[[i]] <- MapQuery(
    anchorset = bal_unt_anchors[[i]], 
    query = list_unt_query[[i]],
    reference = lungs_ref,
    refdata = list(celltype = "celltype"),
    reference.reduction = "pca",
    reduction.model = "umap"
  )
}

bal_unt_lungref <- merge(list_unt_query[[1]], list_unt_query[2:length(list_unt_query)], merge.dr = "ref.umap")

DefaultAssay(bal_unt_lungref) <- "RNA"
bal_unt_lungref <- NormalizeData(bal_unt_lungref)
bal_unt_lungref <- FindVariableFeatures(bal_unt_lungref, selection.method = "vst")
bal_unt_lungref <- ScaleData(bal_unt_lungref)


#-----------------------------------------------------------------------------------------------------------

# SingleR

rData <- read.table("/yerkes-cifs/runs/Analysis/COVID_NHP_Study1/R4_10X_Seuratv4/bulk/new_symbols/rlog_IM_AM_AU.txt", header=T, row.names = NULL)
rownames(rData) <- make.names(rData[,1], unique = TRUE)
rData <- rData[,-1]
head(rData)

pData <- read.table("/yerkes-cifs/runs/Analysis/COVID_NHP_Study1/R4_10X_Seuratv4/bulk/pData.txt", header=T)
pData
bulk <- SummarizedExperiment(assays = list(logcounts=rData), colData = pData)

normdata <- GetAssayData(bal_unt_lungref, assay="RNA", slot = "data")
common <- intersect(row.names(normdata),row.names(bulk))
normdata <- normdata[common,]
bulk <- bulk[common,]

bulk.main <- SingleR(test = normdata, ref = bulk, 
                     labels = bulk$Type, assay.type.ref = "logcounts")

bal_unt_lungref[["Bulk.main"]] <- bulk.main$labels
bal_unt_lungref[["Bulk.main.pruned"]] <- bulk.main$pruned.labels

singleR.markers <- metadata(bulk.main)$de.genes

am.singleR.markers <- unique(unlist(singleR.markers$AM))
non_am.singleR.markers <- unique(unlist(singleR.markers$`non-AM`))

l1 <- sort(unique(unlist(singleR.markers$AM$`non-AM`)))
l2 <- sort(unique(unlist(singleR.markers$`non-AM`$AM)))
write.table(l1,"AM_nonAM_SingelR_markers.txt",quote = F, row.names = F)
write.table(l2,"nonAM_AM_SingelR_markers.txt",quote = F, row.names = F)

#-----------------------------------------------------------------------------------------------------------

# Set up factors

bal_unt_lungref$predicted.celltype <- factor(bal_unt_lungref$predicted.celltype, levels = c("CD163+MRC1+ Mac","CD163+MRC1+TREM2+ Mac","CD163+MRC1- Mac","CD16+ Mono"))


bal_unt_lungref$Time_Group <- paste0(bal_unt_lungref$Timepoint," ",bal_unt_lungref$Group)
bal_unt_lungref$celltype_Group_Time <- paste0(bal_unt_lungref$predicted.celltype,"_",bal_unt_lungref$Group,"_",bal_unt_lungref$Timepoint)
bal_unt_lungref$celltype_Time <- paste0(bal_unt_lungref$predicted.celltype," ",bal_unt_lungref$Timepoint)

bal_unt_lungref$celltype_Time <- factor(bal_unt_lungref$celltype_Time, levels = c("CD163+MRC1- Mac -5dpi","CD163+MRC1- Mac 4dpi",
                                                                                  "CD163+MRC1+TREM2+ Mac -5dpi","CD163+MRC1+TREM2+ Mac 4dpi",
                                                                                  "CD163+MRC1+ Mac -5dpi","CD163+MRC1+ Mac 4dpi",
                                                                                  "CD16+ Mono -5dpi","CD16+ Mono 4dpi"))



saveRDS(bal_unt_lungref,"BAL_untreated_LungRef.rds")
bal_unt_lungref <- readRDS("BAL_untreated_LungRef.rds")

#-----------------------------------------------------------------------------------------------------------

# UMAP

# Fig 4a
DimPlot(bal_unt_lungref, reduction = "ref.umap", group.by =  "Timepoint", repel = TRUE, cols =c("#4d9221","#c51b7d"), pt.size = 1.5, raster = T)  & NoAxes()

# Fig 4b
DimPlot(bal_unt_lungref, reduction = "ref.umap", group.by =  "predicted.celltype", repel = TRUE, split.by = "Timepoint",ncol = 2,
        cols = c("#045a8d","#a50f15","#fb6a4a","#fec44f"),pt.size = 2,raster = T)  & NoAxes()


#-----------------------------------------------------------------------------------------------------------

# Fig 4c

DefaultAssay(bal_unt_lungref) <- "RNA"
Idents(bal_unt_lungref) <- "predicted.celltype"
DotPlot(bal_unt_lungref,features = c("CD163","CD14","FCGR3","MRC1","CHIT1","MARCO","FABP4","TREM2","APOBEC3A","S100A8","LYVE1")) + RotatedAxis()

#-----------------------------------------------------------------------------------------------------------

# Fig 4e

# Proportions

library(tidyr)

# Stacked bar chart - Day
d <- as.data.frame(as.matrix(prop.table(table(bal_unt_lungref$Time_Group,bal_unt_lungref$predicted.celltype),1) * 100))
head(d)
names(d) <- c("Sample","Type","Percentage")
d <- d %>% separate(col=Sample, into = c("Timepoint","Group"), sep=" ", remove = F)
d
ggplot(d, aes(fill=Type, y=Percentage, x=Timepoint)) + 
  geom_bar(position="stack", stat="identity") +
  theme_bw() + ylab("Percentage of Macrophages/Monocytes") + 
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank()) + 
  theme(text = element_text(size=14),
        axis.text.x = element_text(size=12, color = "black"),
        axis.text.y = element_text(size=12, color = "black")) +
  theme(strip.background =element_rect(fill="white")) +
  scale_fill_manual(values= c("#045a8d","#a50f15","#fb6a4a","#fec44f"))

#-----------------------------------------------------------------------------------------------------------

# Fig 4f,g,h

# FeaturePlots

bal_unt_lungref_4dpi <- subset(bal_unt_lungref, cells = row.names(bal_unt_lungref@meta.data[bal_unt_lungref$Timepoint == "4dpi",]))
DefaultAssay(bal_unt_lungref_4dpi) <- "RNA"
bal_unt_lungref_4dpi <- NormalizeData(bal_unt_lungref_4dpi)
bal_unt_lungref_4dpi <- FindVariableFeatures(bal_unt_lungref_4dpi, selection.method = "vst")
bal_unt_lungref_4dpi <- ScaleData(bal_unt_lungref_4dpi)

bal_unt_lungref_4dpi$celltype_Time <- factor(bal_unt_lungref_4dpi$celltype_Time, levels = c("CD163+MRC1- Mac 4dpi","CD163+MRC1+TREM2+ Mac 4dpi","CD163+MRC1+ Mac 4dpi","CD16+ Mono 4dpi"))

saveRDS(bal_unt_lungref_4dpi,"bal_unt_lungref_4dpi.rds")

Idents(bal_unt_lungref_4dpi) <- "celltype_Time"
FeaturePlot(bal_unt_lungref_4dpi, features =c("TNF","IL6","IL10","IFNB1","IL1B"), split.by = "celltype_Time", pt.size = 1.5, cols=c('#DCD9DE',"blue"), max.cutoff = 5, min.cutoff = 0) & NoAxes()
FeaturePlot(bal_unt_lungref_4dpi, features =c("CCL4L1","CXCL10","CXCL3","CXCL8"), split.by = "celltype_Time", pt.size = 1.5, cols=c('#DCD9DE',"blue"), max.cutoff = 5, min.cutoff = 0) & NoAxes()
FeaturePlot(bal_unt_lungref_4dpi, features =c("IFI27","ISG15","ISG20","MX2"), split.by = "celltype_Time", pt.size = 1.5, cols=c('#DCD9DE',"blue"), max.cutoff = 5, min.cutoff = 0) & NoAxes()

#-----------------------------------------------------------------------------------------------------------

# Fig 4i

# contribution CPM

bal_unt_lungref_4dpi_cpm <- NormalizeData(bal_unt_lungref_4dpi, assay = "RNA", normalization.method = "RC", scale.factor = 1e6)
cpm_data <- GetAssayData(bal_unt_lungref_4dpi_cpm, slot = "data", assay = "RNA")
cpm_data_filt <- as.data.frame(cpm_data[c("IL6","IL10","IL1B","IFNB1","TNF","CXCL3","CXCL8","CXCL10","CCL4L1","MX2","ISG15","ISG20","IFI27"),])
write.table(cpm_data_filt,"CPM_contribution.txt", sep = "\t", quote = F)
head(cpm_data_filt,10)[1:10,1:10]
row.names(cpm_data_filt)
total_norm_data_filt <- rowSums(cpm_data_filt)
total_norm_data_filt

am_like <- bal_unt_lungref_4dpi_cpm$predicted.celltype == "CD163+MRC1+ Mac"
am_trem2 <- bal_unt_lungref_4dpi_cpm$predicted.celltype == "CD163+MRC1+TREM2+ Mac"
im_like <- bal_unt_lungref_4dpi_cpm$predicted.celltype == "CD163+MRC1- Mac"
nc <- bal_unt_lungref_4dpi_cpm$predicted.celltype == "CD16+ Mono"

am_like_rowsums <- rowSums(cpm_data_filt[,am_like])
am_trem2_rowsums <- rowSums(cpm_data_filt[,am_trem2])
im_like_rowsums <- rowSums(cpm_data_filt[,im_like])
nc_rowsums <- rowSums(cpm_data_filt[,nc])

all_per <- data.frame("CD163posMRC1pos" = (am_like_rowsums/total_norm_data_filt)*100,
                      "CD163posMRC1posTREM2pos" = (am_trem2_rowsums/total_norm_data_filt)*100,
                      "CD163posMRC1neg" = (im_like_rowsums/total_norm_data_filt)*100, 
                      "CD16pos" = (nc_rowsums/total_norm_data_filt)*100)

all_per$Gene <-  row.names(all_per)

all_per
write.table(all_per,"CPM_percentage.txt",sep="\t",quote = F)

library(reshape2)
d <- melt(all_per)
head(d)
names(d) <- c("Gene","Type1","Percent")
d$Class[d$Gene == "IL6" | d$Gene == "IL10" | d$Gene == "IFNB1" | d$Gene == "TNF" | d$Gene == "IL1B" ] <- "Cytokines"
d$Class[d$Gene == "CXCL10" | d$Gene == "CXCL3" | d$Gene == "CXCL8" | d$Gene == "CCL4L1"] <- "Chemokines"
d$Class[d$Gene == "IFI27" | d$Gene == "ISG15" | d$Gene == "ISG20" | d$Gene == "MX2"] <- "ISG"

d$Type = "NA"
d$Type[d$Type1 == "CD163posMRC1pos"] <- "CD163+MRC1+ Mac"
d$Type[d$Type1 == "CD163posMRC1posTREM2pos"] <- "CD163+MRC1+TREM2+ Mac"
d$Type[d$Type1 == "CD163posMRC1neg"] <- "CD163+MRC1- Mac"
d$Type[d$Type1 == "CD16pos"] <- "CD16+ Mono"

d$Type <- factor(d$Type,levels=c("CD163+MRC1+ Mac","CD163+MRC1+TREM2+ Mac","CD163+MRC1- Mac","CD16+ Mono"))
d$Class <- factor(d$Class,levels=c("Cytokines","Chemokines","ISG"))

ggplot(d,aes(x=Gene,y=Percent,fill=Type)) + 
  geom_bar(position="stack",stat="identity", width = 0.8) +
  facet_wrap(~Class, scales = "free") +
  ylab("Percentage of total normalized\nexpression at 4dpi") +
  theme_bw() + 
  theme(axis.text.x = element_text(color = "black",size = 12, angle = 45, hjust = 1),
        axis.text.y = element_text(color = "black",size = 12),
        axis.title = element_text(color = "black",size = 12),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        strip.background = element_blank(),
        strip.text = element_text(color = "white",size=0),
        panel.border = element_rect(colour = "black", fill = NA)) +
  scale_fill_manual(values=c( c("#045a8d","#a50f15","#fb6a4a","#fec44f"))) 


#-----------------------------------------------------------------------------------------------------------

# Supplementary File 4

# Get DE genes

DefaultAssay(bal_unt_lungref) <- "RNA"
Idents(bal_unt_lungref) <- "Timepoint"
markers_4dpi_base <- FindMarkers(bal_unt_lungref,ident.1 = "4dpi", ident.2 = "-5dpi", test.use = "MAST")
write.table(markers_4dpi_base,"Markers_unt_4dpi_base.txt", quote=F, sep="\t")

DefaultAssay(bal_unt_lungref) <- "RNA"
Idents(bal_unt_lungref) <- "celltype_Time"
markers_trem2_4dpi_base<- FindMarkers(bal_unt_lungref,ident.1 = "CD163+MRC1+TREM2+ Mac 4dpi", ident.2 = "CD163+MRC1+TREM2+ Mac -5dpi", test.use = "MAST")
write.table(markers_trem2_4dpi_base,"Markers_unt_trem2_4dpi_base.txt", quote=F, sep="\t")


DefaultAssay(bal_unt_lungref) <- "RNA"
Idents(bal_unt_lungref) <- "celltype_Time"
markers_am_4dpi_base<- FindMarkers(bal_unt_lungref,ident.1 = "CD163+MRC1+ Mac 4dpi", ident.2 = "CD163+MRC1+ Mac -5dpi", test.use = "MAST")
write.table(markers_am_4dpi_base,"Markers_unt_am_4dpi_base.txt", quote=F, sep="\t")


DefaultAssay(bal_unt_lungref_4dpi) <- "RNA"
Idents(bal_unt_lungref_4dpi) <- "predicted.celltype"
markers_4dpi_trem2_am <- FindMarkers(bal_unt_lungref_4dpi,ident.1 = "CD163+MRC1+TREM2+ Mac", ident.2 = "CD163+MRC1+ Mac", test.use = "MAST")
write.table(markers_4dpi_trem2_am,"Markers_unt_4dpi_trem2_am.txt", quote=F, sep="\t")

DefaultAssay(bal_unt_lungref_4dpi) <- "RNA"
Idents(bal_unt_lungref_4dpi) <- "predicted.celltype"
markers_4dpi_im_am <- FindMarkers(bal_unt_lungref_4dpi,ident.1 = "CD163+MRC1- Mac", ident.2 = "CD163+MRC1+ Mac", test.use = "MAST")
write.table(markers_4dpi_im_am,"Markers_unt_4dpi_im_am.txt", quote=F, sep="\t")

DefaultAssay(bal_unt_lungref_4dpi) <- "RNA"
Idents(bal_unt_lungref_4dpi) <- "predicted.celltype"
markers_4dpi_trem2_im <- FindMarkers(bal_unt_lungref_4dpi,ident.1 = "CD163+MRC1+TREM2+ Mac", ident.2 = "CD163+MRC1- Mac", test.use = "MAST")
write.table(markers_4dpi_trem2_im,"Markers_unt_4dpi_trem2_im.txt", quote=F, sep="\t")



#-----------------------------------------------------------------------------------------------------------

# Fig S6


DefaultAssay(bal_unt_lungref) <- "RNA"

Idents(bal_unt_lungref) <- "celltype_Time"
p1 <- DotPlot(bal_unt_lungref,features=c("CD163","MRC1","CHIT1","FABP4","MARCO","TREM2","ADAMDEC1","S100A8","S100A9","APOBEC3A","GZMB","CXCL3","CXCL8","CCL2","CCL3","CCL5","CCL4L1","CXCL10","IL1B","IL10","TNF","IL6","IFNB1","CCR2","CX3CR1","MKI67","MERTK","FCGR1A","ADGRE1","SIGLEC1","MAMU-DRA","MAMU-DRB1","LYVE1")) + RotatedAxis()

bal_unt_lungref$Bulk_Time <-factor(bal_unt_lungref$Bulk_Time,levels=c("AM -5dpi","AM 4dpi","non-AM -5dpi","non-AM 4dpi"))  
Idents(bal_unt_lungref) <- "Bulk_Time"
p2 <- DotPlot(bal_unt_lungref,features=c("CD163","MRC1","CHIT1","FABP4","MARCO","TREM2","ADAMDEC1","S100A8","S100A9","APOBEC3A","GZMB","CXCL3","CXCL8","CCL2","CCL3","CCL5","CCL4L1","CXCL10","IL1B","IL10","TNF","IL6","IFNB1","CCR2","CX3CR1","MKI67","MERTK","FCGR1A","ADGRE1","SIGLEC1","MAMU-DRA","MAMU-DRB1","LYVE1")) + RotatedAxis()

DefaultAssay(lungs_ref) <- "RNA"
p3 <- DotPlot(lungs_ref,features=c("CD163","MRC1","CHIT1","FABP4","MARCO","TREM2","ADAMDEC1","S100A8","S100A9","APOBEC3A","GZMB","CXCL3","CXCL8","CCL2","CCL3","CCL5","CCL4L1","CXCL10","IL1B","IL10","TNF","IL6","IFNB1","CCR2","CX3CR1","MKI67","MERTK","FCGR1A","ADGRE1","SIGLEC1","MAMU-DRA","MAMU-DRB1","LYVE1")) + RotatedAxis()

p4 <- p3 / p1 / p2
p4
ggsave("FigS6.pdf")

