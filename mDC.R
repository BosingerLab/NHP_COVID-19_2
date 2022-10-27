Idents(int) <- "seurat_clusters"
DotPlot(int, features=c("CCR7","LAMP3","CCL19","MARCKSL1","CD83","IDO1", # Activated
                        "CD1C","FCER1A","CLEC10A", # cDC2
                        "RAB32","CLEC9A","IRF8", # cDC1
                        "MAMU-DRA","MAMU-DRB1")) + coord_flip()

int$mDC = "Others"
int$mDC[int$seurat_clusters == "17"] <- "mDC"
int$mDC[int$seurat_clusters == "22"] <- "Activated mDC"

DefaultAssay(int) <- "RNA"

Idents(int) <- "mDC"
p1 <- DimPlot(int, label = T, cols = c("grey","red","blue"), repel = T,label.size = 4.5, raster = T) & NoAxes() & NoLegend()
p2 <- DotPlot(int, features=c("CCR7","LAMP3","CCL19","MARCKSL1","CD83","IDO1", # Activated
                              "CD1C","FCER1A","CLEC10A", # cDC2
                              "RAB32","CLEC9A","IRF8", # cDC1
                              "MAMU-DRA","MAMU-DRB1")) + coord_flip()


p1 | p2

cytokines <- c("TNF","IL6","IL10","IL1B") 
chemokines <- c("CCL4L1","CXCL10","CXCL3","CXCL8")

int$mDC_group_time <- paste0(int$mDC,"_",int$Group,"_",int$Timepoint)
int$mDC_group_time <- factor(int$mDC_group_time, levels = sort(unique(int$mDC_group_time)))

DefaultAssay(int) <- "RNA"
int$TNFexp <- "TNFneg"
cells_expr <- WhichCells(object = int, expression = (TNF > 0))
int@meta.data[cells_expr,]$TNFexp <- "TNFpos"
table(int$mDC_group_time,int$TNFexp)
prop.table(table(int$mDC_group_time,int$TNFexp),1) * 100

DefaultAssay(int) <- "RNA"
int$IL6exp <- "IL6neg"
cells_expr <- WhichCells(object = int, expression = (IL6 > 0))
int@meta.data[cells_expr,]$IL6exp <- "IL6pos"
table(int$mDC_group_time,int$IL6exp)
prop.table(table(int$mDC_group_time,int$IL6exp),1) * 100


DefaultAssay(int) <- "RNA"
int$IL10exp <- "IL10neg"
cells_expr <- WhichCells(object = int, expression = (IL10 > 0))
int@meta.data[cells_expr,]$IL10exp <- "IL10pos"
table(int$mDC_group_time,int$IL10exp)
prop.table(table(int$mDC_group_time,int$IL10exp),1) * 100


DefaultAssay(int) <- "RNA"
int$IL1Bexp <- "IL1Bneg"
cells_expr <- WhichCells(object = int, expression = (IL1B > 0))
int@meta.data[cells_expr,]$IL1Bexp <- "IL1Bpos"
table(int$mDC_group_time,int$IL1Bexp)
prop.table(table(int$mDC_group_time,int$IL1Bexp),1) * 100
prop.table(table(int$mDC_group_time,int$IL1Bexp),2) * 100

int_untreated <- subset(int, cells = row.names(int@meta.data[int$Group=="Untreated",]))
DefaultAssay(int_untreated) <- "RNA"
int_untreated <- NormalizeData(int_untreated)
int_untreated <- FindVariableFeatures(int_untreated, selection.method = "vst", verbose = F)
int_untreated <- ScaleData(int_untreated, verbose = F)

VlnPlot(int_untreated,cytokines, ncol = 5) 

p3 <- rasterise(VlnPlot(int_untreated,cytokines,stack=T,flip=T) +
            geom_jitter(size=0.01)) + 
  scale_x_discrete(labels = rep(c("-5 dpi","4 dpi"),14)) +
  theme(axis.text.x = element_text(angle = 0, size = 14)) +
  theme(axis.text.y = element_text(angle = 0, size = 14)) +
  theme(axis.text.y = element_text(size=11, color = "black")) +
  geom_vline(xintercept = seq(2.5,6,by=2), linetype="dashed") +
  theme(plot.margin = unit(c(1,1,1,1), "cm")) & NoLegend() 


p4 <- rasterise(VlnPlot(int_untreated,chemokines,stack=T,flip=T) +
                  geom_jitter(size=0.01)) + 
  scale_x_discrete(labels = rep(c("-5 dpi","4 dpi"),14)) +
  theme(axis.text.x = element_text(angle = 0, size = 14)) +
  theme(axis.text.y = element_text(angle = 0, size = 14)) +
  theme(axis.text.y = element_text(size=11, color = "black")) +
  geom_vline(xintercept = seq(2.5,6,by=2), linetype="dashed") +
  theme(plot.margin = unit(c(1,1,1,1), "cm")) & NoLegend() 

d <- as.data.frame(as.matrix(prop.table(table(int$Sample,int$mDC),1) * 100))
head(d)
names(d) <- c("Sample","Type","Percentage")
d <- d %>% separate(col=Sample, into = c("Tissue","Timepoint","Group","AnimalID"), sep="_", remove = F)
d <- d[d$Type != "Others",]
d <- d[d$Group != "Treated",]

ggplot(d, aes(y=Percentage, x=Timepoint)) + 
  geom_line(aes(group=AnimalID), color = "gray") +
  geom_point(size = 4, shape = 21, aes(fill= Type)) +
  theme_bw() +
  ylab("Percentage of mDCs\nout of all BAL cells") + xlab("") +
  scale_y_continuous(expand=expansion(mult=c(0.1,0.2))) +
  scale_fill_manual(values = c("#045a8d","#a50f15","#fb6a4a","#fec44f")) +
  stat_summary(fun= median, fun.min = median, fun.max = median,
               geom = "crossbar", width = 0.75, color = "black", size = 0.5) +
  facet_grid(~Type) +
  theme(strip.background =element_rect(fill="white")) +
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank()) +
  theme(text = element_text(size=14),
        axis.text.x = element_text(color = "black", size = 12),
        axis.text.y = element_text(size=12, color = "black"))


wilcox.test(d[d$Type == "mDC" & d$Timepoint == "-5dpi",]$Percentage,
            d[d$Type == "mDC" & d$Timepoint == "4dpi",]$Percentage,
            paired = T,
            correct = F
)

wilcox.test(d[d$Type == "Activated mDC" & d$Timepoint == "-5dpi",]$Percentage,
            d[d$Type == "Activated mDC" & d$Timepoint == "4dpi",]$Percentage,
            paired = T,
            correct = F
)


int_untreated$mDC_Time <- paste0(int_untreated$mDC,"_",int_untreated$Timepoint)
Idents(int_untreated) <- "mDC_Time"
DefaultAssay(int_untreated) <- ""
mDC_untrt <- FindMarkers(int_untreated, ident.1 = "mDC_4dpi", ident.2 = "mDC_-5dpi", test.use = "MAST")
actmDC_untrt <- FindMarkers(int_untreated, ident.1 = "Activated mDC_4dpi", ident.2 = "Activated mDC_-5dpi", test.use = "MAST")

