library(Seurat)
library(ggplot2)
library(cowplot)
library(ggpubr)
library(dplyr)
library(RColorBrewer)
library(corrplot)
#data <- Read10X_h5("/nfs/scistore18/vicosgrp/melkrewi/Project_snRNA_ovaries_final/2.replicate_1/cellbender/output_FPR_0.3_filtered.h5", use.names = TRUE)
data <- Read10X_h5("/nfs/scistore18/vicosgrp/melkrewi/Project_snRNA_ovaries_with_W/2.replicate_1/cellbender/output_FPR_0.5_filtered.h5", use.names = TRUE)
dataX <- CreateSeuratObject(counts = data)
dataX <- subset(dataX, subset = nFeature_RNA > 10)
total_counts_per_cell <- colSums(dataX@assays$RNA@counts)
mito_genes <- read.csv("mitocondrial_genes.txt")
dataX$percent_mito <- colSums(dataX@assays$RNA@counts[mito_genes$transcript, ])/total_counts_per_cell
feats <- c("nFeature_RNA", "nCount_RNA", "percent_mito")
pdf("QC_dataset_1_round2.pdf")
VlnPlot(dataX, group.by = "orig.ident", features = feats, pt.size = 0.1, ncol = 3) +
    NoLegend()
selected_mito <- WhichCells(dataX, expression = percent_mito < 0.03)
dataX <- subset(dataX, cells = selected_mito)
C <- dataX@assays$RNA@counts
C <- Matrix::t(Matrix::t(C)/Matrix::colSums(C)) * 100
most_expressed <- order(apply(C, 1, median), decreasing = T)[20:1]
boxplot(as.matrix(t(as.matrix(C[most_expressed, ]))), cex = 0.1, las = 1, xlab = "% total count per cell",
    col = (scales::hue_pal())(20)[20:1], horizontal = TRUE)
#dataX <- dataX[!grepl("MSTRG.15825", rownames(dataX)), ]
dataX <- NormalizeData(dataX, normalization.method = "LogNormalize", scale.factor = 10000)
all.genes <- rownames(dataX)
dataX <- ScaleData(dataX,features = all.genes)
dataX <- FindVariableFeatures(dataX, selection.method = "vst", nfeatures = 2000)
dataX <- RunPCA(dataX, features = VariableFeatures(object = dataX))
dataX <- FindNeighbors(dataX, dims = 1:20)
dataX <- FindClusters(dataX, resolution = 0.5)
dataX <- RunUMAP(dataX, dims = 1:20)
DimPlot(dataX, reduction = "umap", label = TRUE, repel = TRUE)
FeaturePlot(dataX, features = c("MSTRG.2397"),order=TRUE)
VlnPlot(dataX, features = c("MSTRG.2397"))
library(dplyr)
library(tibble)
library(Seurat)
export_df <- dataX@meta.data %>%
  rownames_to_column("barcodes") %>%
  select(barcodes, seurat_clusters)
write.csv(export_df, "seurat_clusters_dataset_1_round2_no_DB_removal.txt")
#DimPlot(dataX, reduction = "umap", label = TRUE, repel = TRUE)
#FeaturePlot(dataX, features = c("MSTRG.19052"))
#VlnPlot(dataX, features = c("MSTRG.19052"))
#FeaturePlot(dataX, features = c("MSTRG.19749"))
#VlnPlot(dataX, features = c("MSTRG.19749"))
#FeaturePlot(dataX, features = c("MSTRG.7641"))
#VlnPlot(dataX, features = c("MSTRG.7641"))
#FeaturePlot(dataX, features = c("MSTRG.22546"))
#VlnPlot(dataX, features = c("MSTRG.22546"))
DimPlot(dataX, reduction = "umap", label = TRUE, repel = TRUE)
library(DoubletFinder)
sweep.res.list <- paramSweep_v3(dataX, PCs = 1:20, sct = T)
sweep.stats <- summarizeSweep(sweep.res.list, GT = FALSE)
bcmvn <- find.pK(sweep.stats)
annotations <- dataX@active.ident
homotypic.prop <- modelHomotypic(annotations)
nExp_poi <- round(0.076*length(dataX@active.ident))
nExp_poi.adj <- round(nExp_poi*(1-homotypic.prop))
seu <- doubletFinder_v3(dataX, PCs = 1:10, pN = 0.25, pK = 0.19, nExp = nExp_poi, reuse.pANN = FALSE, sct = FALSE)
DF.name = colnames(seu@meta.data)[grepl("DF.classification", colnames(seu@meta.data))]
cowplot::plot_grid(ncol = 2, DimPlot(seu, group.by = "orig.ident") + NoAxes(),
    DimPlot(seu, group.by = DF.name) + NoAxes())
data.filt = seu[, seu@meta.data[, DF.name] == "Singlet"]
DimPlot(data.filt, reduction = "umap", label = TRUE, repel = TRUE)
library(dplyr)
library(tibble)
library(Seurat)
export_df <- data.filt@meta.data %>%
  rownames_to_column("barcodes") %>%
  select(barcodes, seurat_clusters)
write.csv(export_df, "seurat_clusters_dataset_1_round2.txt")
FeaturePlot(data.filt, features = c("MSTRG.2397"),order=TRUE)
VlnPlot(data.filt, features = c("MSTRG.2397"))
dev.off()
saveRDS(data.filt, "replicate_1.rds")
