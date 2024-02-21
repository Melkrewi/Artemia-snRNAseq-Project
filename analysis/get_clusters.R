library(dplyr)
library(tibble)
library(Seurat)
library(ggplot2)
library(cowplot)
library(ggpubr)
library(dplyr)
library(RColorBrewer)
library(corrplot)
data.filt <- readRDS("data_integrated_harmony_DUBStepR_2.rds")
export_df <- data.filt@meta.data %>%
  rownames_to_column("barcodes") %>%
  select(barcodes, seurat_clusters)
write.csv(export_df, "seurat_clusters_cells_integrated_SAVER.txt")
#export_df <- data.filt@meta.data %>%
#  rownames_to_column("barcodes") %>%
#  select(barcodes, integrated_snn_res.2)
#write.csv(export_df, "seurat_clusters_cells_integrated_standard_workflow_rpca_7_res_2.txt")
#export_df <- data.filt@meta.data %>%
#  rownames_to_column("barcodes") %>%
#  select(barcodes, integrated_snn_res.1.5)
#write.csv(export_df, "seurat_clusters_cells_integrated_standard_workflow_rpca_7_res_1.5.txt")
