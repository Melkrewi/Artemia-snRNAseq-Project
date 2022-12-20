```
library(Seurat)
library(ggplot2)
library(cowplot)
library(ggpubr)
library(dplyr)
library(RColorBrewer)
library(corrplot)
library(DoubletFinder)
data.data <- Read10X_h5("output_filtered.h5", use.names = TRUE)
dataX <- CreateSeuratObject(counts = data.data)
dataX <- subset(dataX, subset = nFeature_RNA > 200)
total_counts_per_cell <- colSums(dataX@assays$RNA@counts)
mito_genes <- read.csv("mitogenes.txt")
dataX$percent_mito <- colSums(dataX@assays$RNA@counts[mito_genes$transcript, ])/total_counts_per_cell
feats <- c("nFeature_RNA", "nCount_RNA", "percent_mito")
pdf("percent_mito.pdf")
VlnPlot(dataX, group.by = "orig.ident", features = feats, pt.size = 0.1, ncol = 3) +
    NoLegend()
dev.off()
selected_mito <- WhichCells(dataX, expression = percent_mito < 0.05)
dataX <- subset(dataX, cells = selected_mito)
C <- dataX@assays$RNA@counts
C <- Matrix::t(Matrix::t(C)/Matrix::colSums(C)) * 100
most_expressed <- order(apply(C, 1, median), decreasing = T)[20:1]
pdf("most_expressed_genes.pdf")
boxplot(as.matrix(t(as.matrix(C[most_expressed, ]))), cex = 0.1, las = 1, xlab = "% total count per cell",
    col = (scales::hue_pal())(20)[20:1], horizontal = TRUE)
dev.off()
dataX <- dataX[!grepl("MSTRG.15825", rownames(dataX)), ]
dataX <- NormalizeData(dataX, normalization.method = "LogNormalize", scale.factor = 10000)
all.genes <- rownames(dataX)
dataX <- ScaleData(dataX,features = all.genes)
dataX <- FindVariableFeatures(dataX, selection.method = "vst", nfeatures = 2000)
dataX <- RunPCA(dataX, features = VariableFeatures(object = dataX))
dataX <- FindNeighbors(dataX, dims = 1:20)
dataX <- FindClusters(dataX, resolution = 0.5)
dataX <- RunUMAP(dataX, dims = 1:20)
pdf("umap_before_doublet_removal.pdf")
DimPlot(dataX, reduction = "umap", label = TRUE, repel = TRUE)
dev.off()
sweep.res.list <- paramSweep_v3(dataX, PCs = 1:20, sct = T)
sweep.stats <- summarizeSweep(sweep.res.list, GT = FALSE)
bcmvn <- find.pK(sweep.stats)
annotations <- dataX@active.ident
homotypic.prop <- modelHomotypic(annotations)
nExp_poi <- round(0.076*length(dataX@active.ident))
nExp_poi.adj <- round(nExp_poi*(1-homotypic.prop))
seu <- doubletFinder_v3(dataX, PCs = 1:10, pN = 0.25, pK = 0.19, nExp = nExp_poi, reuse.pANN = FALSE, sct = FALSE)
DF.name = colnames(seu@meta.data)[grepl("DF.classification", colnames(seu@meta.data))]
pdf("orig_vs_doublets.pdf")
cowplot::plot_grid(ncol = 2, DimPlot(seu, group.by = "orig.ident") + NoAxes(),
    DimPlot(seu, group.by = DF.name) + NoAxes())
dev.off()
data.filt = seu[, seu@meta.data[, DF.name] == "Singlet"]
pdf("umap_after_doublet_removal.pdf")
DimPlot(data.filt, reduction = "umap", label = TRUE, repel = TRUE)
dev.off()
data.filt.markers <- FindAllMarkers(data.filt, only.pos = TRUE, min.pct = 0.25, logfc.threshold = 0.25)
write.table(data.filt.markers,"markers_all_7082.txt")
data.filt.markers %>%
    group_by(cluster) %>%
    top_n(n = 50, wt = avg_log2FC) -> top50
write.table(top50,"markers_top50_cellbender_22_11_2022_06_02_7082.txt")
data.filt.markers %>%
    group_by(cluster) %>%
    top_n(n = 100, wt = avg_log2FC) -> top100
write.table(top100,"markers_top100_cellbender_22_11_2022_06_02_7082.txt")
library(dplyr)
library(tibble)
library(Seurat)
export_df <- data.filt@meta.data %>%
  rownames_to_column("barcodes") %>%
  select(barcodes, seurat_clusters)
write.csv(export_df, "seurat_clusters_7082_Cells.txt")

data.filt <- readRDS("data_filt_small.rds")
savehistory(file = ".Rhistory_22_11_2022_7082")
write.table(data.filt.markers,"markers_all_7082.txt")
data.filt.markers %>%
    group_by(cluster) %>%
    top_n(n = 50, wt = avg_log2FC) -> top50
write.table(top50,"markers_top50_cellbender_22_11_2022_06_02_7082.txt")
data.filt.markers %>%
    group_by(cluster) %>%
    top_n(n = 100, wt = avg_log2FC) -> top100
write.table(top100,"markers_top100_cellbender_22_11_2022_06_02_7082.txt")
library(dplyr)
library(tibble)
library(Seurat)
export_df <- data.filt@meta.data %>%
  rownames_to_column("barcodes") %>%
  select(barcodes, seurat_clusters)
write.csv(export_df, "seurat_clusters_7082_Cells.txt")
saveRDS(data.filt, "data_filt_small.rds")

savehistory(file = ".Rhistory_22_11_2022_7082")
```
