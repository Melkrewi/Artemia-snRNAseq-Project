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
new.cluster.ids <- c("Follicle cells","Germ cells A","Ovarian muscle cells","Escort cells","Tracheal cells","Prefollicle cells","Germ cells B")
names(new.cluster.ids) <- levels(seurat_obj)
seurat_obj <- RenameIdents(seurat_obj, new.cluster.ids)
seurat_obj$seurat_clusters <- Idents(seurat_obj)
levels(data) <- sort(levels(seurat_obj))
MEs <- GetMEs(seurat_obj, wgcna_name='modules')
modules <- GetModules(seurat_obj)
mods <- levels(modules$module)
mods <- mods[mods!='grey']
meta <- seurat_obj@meta.data
seurat_obj@meta.data <- cbind(meta, MEs)
p <- DotPlot(seurat_obj, features=c('turquoise','green','brown','black','tan'), group.by = 'seurat_clusters')
# flip the x/y axes, rotate the axis labels, and change color scheme:
p <- p + RotatedAxis() + scale_color_gradient2(high='red', mid='grey95', low='blue')
# plot output
pdf("modules_germline_specific.pdf",width=6,height=4)
p
dev.off()
