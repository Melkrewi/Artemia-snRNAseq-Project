library(leidenbase)
library(Seurat)
library(ggplot2)
library(cowplot)
library(ggpubr)
library(dplyr)
library(RColorBrewer)
library(corrplot)
library(SeuratWrappers)
library(monocle3)
library(Seurat)
# plotting and data science packages
library(tidyverse)
library(cowplot)
library(patchwork)
# co-expression network analysis packages:
library(WGCNA)
library(hdWGCNA)
# using the cowplot theme for ggplot
theme_set(theme_cowplot())
# set random seed for reproducibility
set.seed(12345)
seurat_obj <- readRDS('modules_all_celltypes.rds')
#seurat_obj <- seurat_obj[,  seurat_obj$seurat_clusters %in% c("1","6")]
cds <- as.cell_data_set(seurat_obj)

# run the monocle clustering
cds <- cluster_cells(cds, reduction_method='UMAP')

# learn graph for pseudotime
cds <- learn_graph(cds)

p1 <- plot_cells(
  cds = cds,
  color_cells_by = "seurat_clusters",
  show_trajectory_graph = TRUE,
  label_principal_points = TRUE
) 

# plot the UMAP partitions from the clustering algorithm
p2 <- plot_cells(
  cds = cds,
  color_cells_by = "partition",
  show_trajectory_graph = FALSE
)

p1+p2
principal_node <- 'Y_2'
cds <- order_cells(cds,root_pr_nodes = principal_node)

# add pseudotime to seurat object:
seurat_obj$pseudotime <- pseudotime(cds)
seurat_obj$germline_pseudotime <- ifelse(seurat_obj$seurat_clusters %in% c("1", "6"), seurat_obj$pseudotime, NA)

library(viridis)

seurat_obj$UMAP1 <- seurat_obj@reductions$umap@cell.embeddings[,1]
seurat_obj$UMAP2 <- seurat_obj@reductions$umap@cell.embeddings[,2]

p1 <- seurat_obj@meta.data %>%
  ggplot(aes(x=UMAP1, y=UMAP2, color=germline_pseudotime)) +
  ggrastr::rasterise(geom_point(size=1), dpi=500, scale=0.75) +
  coord_equal() +
  scale_color_gradientn(colors=viridis(256), na.value='grey') +
  umap_theme()
p1
pdf("pseudo_time_whole_data_2.pdf",width=12,height=8)
p2 <- plot_cells(
  cds = cds,
  color_cells_by = "seurat_clusters",
  show_trajectory_graph = TRUE,
  label_principal_points = TRUE) 
p2 + p1
dev.off()
saveRDS(seurat_obj, file='germline_pseudotime.rds')
