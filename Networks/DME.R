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
group3 <- seurat_obj@meta.data %>% subset(seurat_clusters %in% c(1,6)) %>% rownames
group4 <- seurat_obj@meta.data %>% subset(seurat_clusters %in% c(0,2,4,5,3)) %>% rownames
DMEs_g_vs_s <- FindDMEs(
  seurat_obj,
  barcodes1 = group3,
  barcodes2 = group4,
  test.use='wilcox',
  wgcna_name='modules'
)
write.csv(DMEs_g_vs_s, row.names=FALSE, quote=FALSE, file='DMEs_g_vs_s_all.csv')
