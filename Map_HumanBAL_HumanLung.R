library("Seurat")
library("ggplot2")
library("dplyr")
library("stringr")
library("SingleR")

ref <- BlueprintEncodeData()

setwd("/data/Analysis_BAL_Lung_Human_2022")

options(future.globals.maxSize= 33554432000) # 32000*1024^2

list_bal_human <- readRDS("list_bal_human_macromono_DietSeurat.rds")

human_ref <- readRDS("lung_int_kropski_VUH.rds")
DefaultAssay(human_ref) <- "integrated"
human_ref <- RunUMAP(human_ref, dims = 1:15, reduction = "pca", return.model = TRUE)

list_query_human <- list_bal_human
DefaultAssay(human_ref) <- "integrated"

bal_human_anchors <- list()
for (i in 1:length(list_query_human)) {
  bal_human_anchors[[i]] <- FindTransferAnchors(
    reference = human_ref,
    query = list_query_human[[i]],
    reference.reduction = "pca",
    dims = 1:15
  )
}

for (i in 1:length(list_query_human)) {
  list_query_human[[i]] <- MapQuery(
    anchorset = bal_human_anchors[[i]], 
    query = list_query_human[[i]],
    reference = human_ref,
    refdata = list(celltype = "Type"),
    reference.reduction = "pca",
    reduction.model = "umap"
  )
}


human_lungref <- merge(list_query_human[[1]], list_query_human[2:length(list_query_human)], merge.dr = "ref.umap")

DefaultAssay(human_lungref) <- "RNA"
human_lungref <- NormalizeData(human_lungref)
human_lungref <- FindVariableFeatures(human_lungref, selection.method = "vst")
human_lungref <- ScaleData(human_lungref)


human_lungref$predicted.celltype <- factor(human_lungref$predicted.celltype, levels = c("FABP4hi","SPP1hi","FCN1hi","Proliferating"))
Idents(human_lungref) <- "Type"

DimPlot(human_lungref, group.by = "predicted.celltype", split.by = "Group",
        cols = c("#74a9cf","#d7301f","#fcbba1","#238b45"),pt.size = 2, raster=T) & NoAxes()

Idents(human)
FeaturePlot(human_lungref, features = c("FABP4hi","SPP1hi","FCN1hi","Proliferating"),  
            reduction = "ref.umap", cols = c("lightgrey", "darkred"), ncol = 2) & theme(plot.title = element_text(size = 10)) & NoAxes()

saveRDS(human_lungref,"human_lungref.rds")

####################################################################################################################

human_lungref <- readRDS("human_lungref.rds")

human_lungref$Group_Type <- paste0(human_lungref$Group,"\n",human_lungref$predicted.celltype)
Idents(human_lungref) <- "Group_Type"

human_lungref$Group_Type <- factor(human_lungref$Group_Type, levels = c("Healthy\nFABP4hi","Healthy\nSPP1hi","Healthy\nFCN1hi","Healthy\nProliferating",
                                                                        "Mild\nFABP4hi","Mild\nSPP1hi","Mild\nFCN1hi","Mild\nProliferating",
                                                                        "Severe\nFABP4hi","Severe\nSPP1hi","Severe\nFCN1hi","Severe\nProliferating"))
levels(human_lungref$Group_Type)
Idents(human_lungref) <- "Group_Type"
p1 <- DotPlot(human_lungref,features = c("TNF","IL6","IL10","IFNB1","IL1B"), assay = "RNA") + coord_flip()
p2 <- DotPlot(human_lungref,features = c("CCL4L1","CXCL10","CXCL3","CXCL8","CCR2"), assay = "RNA") + coord_flip()
p3 <- DotPlot(human_lungref,features = c("IFI27","ISG15","ISG20","MX2"), assay = "RNA") + coord_flip()

p1/p2/p3


Idents(human_lungref) <- "predicted.celltype"
genes2 <- c("CD163","MRC1","MARCO","FABP4","INHBA","TREM2","SPP1","MERTK","LGMN","SIGLEC10","IL1B","FCN1","VCAN","S100A8","S100A9")
DotPlot(human_lungref, features = genes2, assay = "RNA") + coord_flip()


# contribution CPM

get_contribution <- function(obj,type){
  obj <- NormalizeData(obj, assay = "RNA", normalization.method = "RC", scale.factor = 1e6, verbose = F)
  cpm_data <- GetAssayData(obj, slot = "data", assay = "RNA")
  cpm_data_filt <- as.data.frame(cpm_data[c("IL6","IL10","IL1B","IFNB1","TNF","CXCL3","CXCL8","CXCL10","CCR2","MX2","ISG15","ISG20","IFI27"),])
  write.table(cpm_data_filt,paste0("CPM_contribution_",type,".txt"), sep = "\t", quote = F)
  
  row.names(cpm_data_filt)
  total_norm_data_filt <- rowSums(cpm_data_filt)
  
  fabp4 <- obj$predicted.celltype == "FABP4hi"
  spp1 <- obj$predicted.celltype == "SPP1hi"
  fcn1 <- obj$predicted.celltype == "FCN1hi"
  prolif <- obj$predicted.celltype == "Proliferating"
  
  fabp4_rowsums <- rowSums(cpm_data_filt[,fabp4])
  spp1_rowsums <- rowSums(cpm_data_filt[,spp1])
  fcn1_rowsums <- rowSums(cpm_data_filt[,fcn1])
  prolif_rowsums <- rowSums(cpm_data_filt[,prolif])
  
  all_per <- data.frame("FABP4hi" = (fabp4_rowsums/total_norm_data_filt)*100,
                        "SPP1hi" = (spp1_rowsums/total_norm_data_filt)*100,
                        "FCN1hi" = (fcn1_rowsums/total_norm_data_filt)*100,
                        "Proliferating" = (prolif_rowsums/total_norm_data_filt)*100)
  
  all_per$Gene <-  row.names(all_per)
  write.table(all_per,paste0("CPM_percentage_",type,".txt"),sep="\t",quote = F)
  
  library(reshape2)
  d <- melt(all_per)
  head(d)
  names(d) <- c("Gene","Type1","Percent")
  d$Class[d$Gene == "IL6" | d$Gene == "IL10" | d$Gene == "IFNB1" | d$Gene == "TNF" | d$Gene == "IL1B" ] <- "Cytokines"
  d$Class[d$Gene == "CXCL10" | d$Gene == "CXCL3" | d$Gene == "CXCL8" | d$Gene == "CCR2"] <- "Chemokines"
  d$Class[d$Gene == "IFI27" | d$Gene == "ISG15" | d$Gene == "ISG20" | d$Gene == "MX2"] <- "ISG"
  
  d$Class <- factor(d$Class,levels=c("Cytokines","Chemokines","ISG"))
  
  ggplot(d,aes(x=Gene,y=Percent,fill=Type1)) + 
    geom_bar(position="stack",stat="identity", width = 0.8) +
    facet_wrap(~Class, scales = "free") +
    ylab(paste0("Percentage of total normalized\nexpression in ",type," samples")) +
    theme_bw() + 
    theme(axis.text.x = element_text(color = "black",size = 12, angle = 45, hjust = 1),
          axis.text.y = element_text(color = "black",size = 12),
          axis.title = element_text(color = "black",size = 12),
          panel.grid.major = element_blank(),
          panel.grid.minor = element_blank(),
          strip.background = element_blank(),
          strip.text = element_text(color = "white",size=0),
          panel.border = element_rect(colour = "black", fill = NA)) +
    scale_fill_manual(values=c( c("#74a9cf","#d7301f","#fcbba1","#238b45"))) + ggtitle(type)
}

bal_control_cpm <- subset(human_lungref, cells = row.names(human_lungref@meta.data[human_lungref$Group == "Healthy",]))
bal_mild_cpm <- subset(human_lungref, cells = row.names(human_lungref@meta.data[human_lungref$Group == "Mild",]))
bal_severe_cpm <- subset(human_lungref, cells = row.names(human_lungref@meta.data[human_lungref$Group == "Severe",]))

p1 <- get_contribution(bal_control_cpm,"Healthy")
p2 <- get_contribution(bal_mild_cpm,"Mild")
p3 <- get_contribution(bal_severe_cpm,"Severe")

p1/p2/p3

# Proportions

d <- as.data.frame(prop.table(table(human_lungref$Sample,human_lungref$predicted.celltype),1) * 100)
names(d) <- c("Sample","CellType","Percent")
d
d$Group <- "Severe"
d$Group[grepl("C51|C52|C100",d$Sample)] <- "Healthy"
d$Group[grepl("C141|C142|C144",d$Sample)] <- "Mild"


df <- data.frame()
for(x in unique(d$CellType)){
  print(x)
  mild <- d[d$CellType == x & d$Group == "Mild",]$Percent
  severe <- d[d$CellType == x & d$Group == "Severe",]$Percent
  healthy <- d[d$CellType == x & d$Group == "Healthy",]$Percent
  
  df[x,"Healthy vs Mild"] <- wilcox.test(healthy,mild, paired = F, correct = F)$p.value
  df[x,"Healthy vs Severe"] <- wilcox.test(healthy,severe, paired = F, correct = F)$p.value
  df[x,"Severe vs Mild"] <- wilcox.test(severe,mild, paired = F, correct = F)$p.value
  }
df

ggplot(d,aes(x=CellType,y=Percent, group = Group)) +
  geom_jitter(size = 3, shape = 21,aes(group=Group,fill=Group), position = position_dodge(width = 1)) +
  scale_fill_manual(values = c("#fdd0a2","#f16913","#a63603")) +
  stat_summary(fun= median, geom = "crossbar", width = 0.75, color = "black", size = 0.5, position = position_dodge(width = 1)) +
  ylab("Percent of macrophages/monocytes") +
  theme_bw() +
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank()) + 
  theme(text = element_text(size=14),
        axis.text.x = element_text(size=12, color = "black"),
        axis.text.y = element_text(size=12, color = "black")) 

  
