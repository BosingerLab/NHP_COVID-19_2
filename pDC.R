bal_int_rhesus <- readRDS("/data/BAL/BAL_int_filt_Myeloid.rds")

DefaultAssay(bal_int_rhesus) <- "RNA"

p1 <- DimPlot(bal_int_rhesus, group.by = "seurat_clusters", label = T, raster = T)
p2 <- DimPlot(bal_int_rhesus, group.by = "BP.main", label = T, raster = T)

Idents(bal_int_rhesus) <- "seurat_clusters"
p3 <- DotPlot(bal_int_rhesus, assay = "RNA",features = c("CD19","MS4A1","CD3D","CD3E","CD4","LILRA4","IL3RA","CLEC4C","FLT3","CCDC50","GZMB","JCHAIN","TLR7","TLR9","IRF7","IRF8","TIGIT",
                                           "IFNA13","IFNA8","IFNA6","IFNA14","IFNA10","IFNB1","MAMU-DRA","MAMU-DRB1","MARCO","MRC1","CD163")) + coord_flip()

(p1 + p2 & NoLegend() & NoAxes())/p3 + plot_layout(heights = c(1,2)) 

# Proportion

library(tidyr)
d <- as.data.frame(prop.table(table(bal_int_rhesus$Sample,bal_int_rhesus$seurat_clusters),1) * 100)
d
d <-d[d$Var2 ==11,]
d
names(d) <- c("Sample","cluster","Percentage")
d <- d %>% separate(col="Sample", into = c("Tissue","Timepoint","Group","AnimalID"), sep = "_", remove = FALSE, )
d$Sample <- gsub("BAL_","",d$Sample)
d$Sample <- gsub("_"," ",d$Sample)
d$Group[d$Group == "Treated"] <- "+Baricitinib"
d$Group <- factor(d$Group, levels = c("Untreated","+Baricitinib"))

p1 <- ggplot(d,aes(x=Timepoint,y=Percentage)) +
  geom_line(aes(group = AnimalID),position = position_dodge(width = 0.2)) +
  geom_jitter(aes(fill=AnimalID), shape = 21, size = 4,position = position_dodge(width = 0.2)) + 
  facet_wrap(~Group) +
  theme_bw() +
  ylab("Percentage of pDC\nout of all cell types") +
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank()) + 
  theme(text = element_text(size=14),
        axis.text.x = element_text(size=12, color = "black"),
        axis.text.y = element_text(size=12, color = "black"),
        strip.background =element_rect(fill="white"))   
  

# Absolute counts

d <- as.data.frame(table(bal_int_rhesus$Sample,bal_int_rhesus$seurat_clusters))
d <-d[d$Var2 ==11,]
d

names(d) <- c("Sample","cluster","Percentage")
d <- d %>% separate(col="Sample", into = c("Tissue","Timepoint","Group","AnimalID"), sep = "_", remove = FALSE, )
d$Sample <- gsub("BAL_","",d$Sample)
d$Sample <- gsub("_"," ",d$Sample)
d$Group[d$Group == "Treated"] <- "+Baricitinib"
d$Group <- factor(d$Group, levels = c("Untreated","+Baricitinib"))

p2 <- ggplot(d,aes(x=Timepoint,y=Percentage)) +
  geom_line(aes(group = AnimalID),position = position_dodge(width = 0.2)) +
  geom_jitter(aes(fill=AnimalID), shape = 21, size = 4,position = position_dodge(width = 0.2)) + 
  facet_wrap(~Group) +
  theme_bw() +
  ylab("Absolute number of pDC") +
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank()) + 
  theme(text = element_text(size=14),
        axis.text.x = element_text(size=12, color = "black"),
        axis.text.y = element_text(size=12, color = "black"),
        strip.background =element_rect(fill="white"))   

p2 + p1


