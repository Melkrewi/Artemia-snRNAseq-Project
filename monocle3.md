# monocle3 pseudotime 
```
library(Seurat)
library(ggplot2)
library(cowplot)
library(ggpubr)
library(dplyr)
library(RColorBrewer)
library(corrplot)
data.filt <- readRDS("data_filt_small.rds")
library(SeuratWrappers)
library(monocle3)
germline_cells <- data.filt[,  data.filt$seurat_clusters %in% c("2","3","7","10")]
cds <- as.cell_data_set(germline_cells)
recreate.partitions <- c(rep(1, length(cds@colData@rownames)))
names(recreate.partitions) <- cds@colData@rownames
recreate.partitions <- as.factor(recreate.partitions)
recreate.partitions
cds@clusters@listData[["UMAP"]][["partitions"]] <- recreate.partitions
germline_cells@active.ident
list.cluster <- germline_cells@active.ident
cds@clusters@listData[["UMAP"]][["clusters"]] <- list.cluster
cds@int_colData@listData[["reducedDims"]]@listData[["UMAP"]] <- germline_cells@reductions$umap@cell.embeddings
cluster.before.traj <-plot_cells(cds, color_cells_by = "cluster", label_groups_by_cluster = F, 
           group_label_size = 5) + theme(legend.position = "right")
cluster.before.traj
cds <- learn_graph(cds, use_partition = F)

plot_cells(cds, color_cells_by = "cluster", label_groups_by_cluster = F,
           label_branch_points = T, label_roots = T, label_leaves = F,
           group_label_size = 5)
cds <- order_cells(cds, reduction_method = "UMAP", root_cells = colnames(cds[, clusters(cds) == 10]))
pdf("pseudotime_10_start.pdf")
plot_cells(cds, color_cells_by = "pseudotime", label_groups_by_cluster = T,
           label_branch_points = T, label_roots = F, label_leaves = F)
dev.off()
```
