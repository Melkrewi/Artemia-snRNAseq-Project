library(Seurat) 
library(ggpubr) 
library(dplyr)
library(ggplot2) 
library(cowplot) 
library(RColorBrewer) 
library(corrplot) 
library(data.table)
set.seed(1234)
data <- readRDS("data_integrated_harmony_DUBStepR_2.rds")
new.cluster.ids <- c("Early Germline cells","Ovarian muscle cells","Follicle cells A/Escort cells","Tracheal cells/Follicle cells C","Follicle cells B","Prefollicle cells","Unannotated","Late Germline cells")
names(new.cluster.ids) <- levels(data)
data <- RenameIdents(data, new.cluster.ids)
data$seurat_clusters <- Idents(data)
levels(data) <- sort(levels(data)) 
#png("umap_renamed_2.png",width=800,height=800)
pdf("umap_renamed_2.pdf",width=10,height=7)
DimPlot(data, reduction = "umap", label = TRUE, repel = TRUE)
dev.off()




