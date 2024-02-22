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
seurat_obj <- readRDS('../../../data_integrated_harmony_DUBStepR_3_2.rds')
seurat_obj <- SetupForWGCNA(
  seurat_obj,
  gene_select = "fraction",
  fraction = 0.05,
  wgcna_name = 'modules'
)

# construct metacells 
seurat_obj <- MetacellsByGroups(
  seurat_obj = seurat_obj,
  group.by = c("seurat_clusters","replicates"), # specify the columns in seurat_obj@meta.data to group by
  reduction = 'umap', # select the dimensionality reduction to perform KNN on
  k = 25, # nearest-neighbors parameter
  max_shared = 10, # maximum number of shared cells between two metacells
  ident.group = 'seurat_clusters' # set the Idents of the metacell seurat object
)
seurat_obj <- NormalizeMetacells(seurat_obj)

# setup expression matrix
seurat_obj <- SetDatExpr(
  seurat_obj,
  group.by='seurat_clusters',
  group_name = c('0','1','2','3','4','5','6'),assay = 'RNA', # using RNA assay
  slot = 'data' # using normalized data
)
seurat_obj <- TestSoftPowers(seurat_obj)
seurat_obj <- ConstructNetwork(
    seurat_obj, 
    tom_name='modules', 
    overwrite_tom=TRUE
)

# compute module eigengenes & connectivity
seurat_obj <- ModuleEigengenes(seurat_obj)
seurat_obj <- ModuleConnectivity(seurat_obj)
saveRDS(seurat_obj, file='modules_all_celltypes.rds')
# plot dendro
PlotDendrogram(seurat_obj, main='hdWGCNA Dendrogram')
MEs <- GetMEs(seurat_obj)
modules <- GetModules(seurat_obj)
mods <- levels(modules$module)
mods <- mods[mods!='grey']
meta <- seurat_obj@meta.data
seurat_obj@meta.data <- cbind(meta, MEs)
# make dotplot
p <- DotPlot(
  seurat_obj,
  group.by='seurat_clusters',
  features = rev(mods)
) + RotatedAxis() +
  scale_color_gradient2(high='red', mid='grey95', low='blue') + xlab('') + ylab('') +
  theme(
    plot.title = element_text(hjust = 0.5),
    axis.line.x = element_blank(),
    axis.line.y = element_blank(),
    panel.border = element_rect(colour = "black", fill=NA, size=1)
  )

p
write.table(modules, "modules_seurat_cluster_all.txt",sep="\t")
