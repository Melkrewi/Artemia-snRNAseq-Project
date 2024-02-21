library(harmony)
library(Seurat) 
library(ggpubr) 
library(dplyr)
library(ggplot2) 
library(cowplot) 
library(RColorBrewer) 
library(corrplot) 
library(data.table)
library(Signac)
library(scCustomize)
library(DUBStepR)
set.seed(1234)
data1 <- Read10X_h5("/nfs/scistore18/vicosgrp/melkrewi/Project_snRNA_ovaries_with_W/2.replicate_1/cellbender/output_FPR_0.5_filtered.h5", use.names = TRUE)
barcodes1 <- as.character(read.csv("../seurat_clusters_dataset_1_round2.txt", header = TRUE)$barcodes)
data1 <- data1[, which((colnames(data1) %in% barcodes1)==TRUE)]
data2 <- Read10X_h5("/nfs/scistore18/vicosgrp/melkrewi/Project_snRNA_ovaries_with_W/3.replicate_2/cellbender/output_FPR_0.5_filtered.h5", use.names = TRUE)
barcodes2 <- as.character(read.csv("../seurat_clusters_dataset_2_round2.txt", header = TRUE)$barcodes)
data2 <- data2[, which((colnames(data2) %in% barcodes2)==TRUE)]
data3 <- Read10X_h5("/nfs/scistore18/vicosgrp/melkrewi/Project_snRNA_ovaries_with_W/4.replicate_3/cellbender/output_FPR_0.5_filtered.h5", use.names = TRUE)
barcodes3 <- as.character(read.csv("../seurat_clusters_dataset_3_round2.txt", header = TRUE)$barcodes)
data3 <- data3[, which((colnames(data3) %in% barcodes3)==TRUE)]
#data4 <- Read10X_h5("/nfs/scistore18/vicosgrp/melkrewi/Project_snRNA_ovaries_with_W/5.replicate_4/cellbender/output_FPR_0.4_filtered.h5", use.names = TRUE)
data4 <- Read10X_h5("/nfs/scistore18/vicosgrp/melkrewi/Project_snRNA_ovaries_with_W/5.replicate_4/cellbender/round2/round3/round4/output_FPR_0.25_filtered.h5", use.names = TRUE)
barcodes4 <- as.character(read.csv("seurat_clusters_dataset_4_round2.txt", header = TRUE)$barcodes)
data4 <- data4[, which((colnames(data4) %in% barcodes4)==TRUE)]
data <- CreateSeuratObject(counts = cbind(data1, data2, data3, data4), project = "Artemia", min.cells = 5) 
data@meta.data$replicates <- c(rep("replicate_1", ncol(data1)), rep("replicate_2", ncol(data2)),rep("replicate_3", ncol(data3)),rep("replicate_4", ncol(data4)))
data <- subset(data, subset = nFeature_RNA > 200 & nFeature_RNA < 20000)
pdf("results_with_harmony_DUBStepR_2.pdf",width=10,height=7)
data <- Seurat::NormalizeData(data,verbose = FALSE)
dubstepR.out <- DUBStepR(input.data = data@assays$RNA@data, min.cells = 0.05*ncol(data), optimise.features = T, k = 10, num.pcs = 20, error = 0,species = "Artemia")
data@assays$RNA@var.features <- dubstepR.out$optimal.feature.genes
data<- ScaleData(data,verbose = FALSE)  %>% RunPCA(pc.genes = data@var.genes, verbose = FALSE,dims = 1:20)
options(repr.plot.height = 5, repr.plot.width = 12)
p1 <- DimPlot(object = data, reduction = "pca", pt.size = .1, group.by = "replicates")
p2 <- VlnPlot(object = data, features = "PC_1", group.by = "replicates", pt.size = .1)
options(repr.plot.height = 2.5, repr.plot.width = 6)
data <- data %>% 
    RunHarmony("replicates", plot_convergence = TRUE,theta = 0.05)
harmony_embeddings <- Embeddings(data, 'harmony')
harmony_embeddings[1:5, 1:5]
options(repr.plot.height = 5, repr.plot.width = 12)
p1 <- DimPlot(object = data, reduction = "harmony", pt.size = .1, group.by = "replicates")
p2 <- VlnPlot(object = data, features = "harmony_1", group.by = "replicates", pt.size = .1)
plot_grid(p1,p2)
data <- data %>% 
    RunUMAP(reduction = "harmony",dims = 1:20, seed.use = 2023) %>% 
    FindNeighbors(reduction = "harmony", dims = 1:20) %>% 
    FindClusters(resolution = c(2,1.5,0.5,0.2,0.1)) %>% 
    identity()
#min.dist=0.2
pdf("results_with_all_stable_harmony_DUBStepR_2.pdf",width=10,height=7)
DimPlot(data, group.by = "replicates")
DefaultAssay(data) <- "RNA"
FeaturePlot(data, features = c("MSTRG.2397"),order=TRUE,split.by = "replicates")
FeaturePlot(data, features = c("MSTRG.17673"),order=TRUE,split.by = "replicates")
VlnPlot(data, features = c("MSTRG.2397"),split.by = "replicates")
VlnPlot(data, features = c("MSTRG.17673"),split.by = "replicates")
DimPlot(data, reduction = "umap", label = TRUE, label.size = 2.5, repel = TRUE) + ggtitle("RNA")
dev.off()
saveRDS(data, "data_integrated_harmony_DUBStepR_2.rds")
