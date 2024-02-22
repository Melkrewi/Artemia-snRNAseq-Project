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
seurat_obj <- readRDS('germline_pseudotime.rds')
p  <- PlotModuleTrajectory(
    seurat_obj,
    pseudotime_col = 'germline_pseudotime',wgcna='modules')
pdf("module_dynamics.pdf")
p
dev.off()
