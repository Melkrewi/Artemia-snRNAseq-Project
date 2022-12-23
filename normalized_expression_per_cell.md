# get csv of normalized expression (genes x barcodes) 
I used R/4.2.1
```
library(Seurat)
library(ggplot2)
library(cowplot)
library(ggpubr)
library(dplyr)
library(RColorBrewer)
library(corrplot)
data.filt <- readRDS("data_filt_small.rds")
library(scrattch.io)
write_dgCMatrix_csv(data.filt@assays$RNA@data, "normalized_expression.csv")
```
