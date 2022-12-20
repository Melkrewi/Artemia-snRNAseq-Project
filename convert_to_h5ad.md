```
library(Seurat)
library(ggplot2)
library(cowplot)
library(ggpubr)
library(dplyr)
library(RColorBrewer)
library(corrplot)
data.filt <- readRDS("data_filt_small.rds")
library(SeuratDisk)
SaveH5Seurat(data.filt, filename = "data_filt.h5Seurat")
Convert("data_filt.h5Seurat", dest = "h5ad")
```
