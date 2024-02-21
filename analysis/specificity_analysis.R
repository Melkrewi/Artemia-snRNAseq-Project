library(scMiko)
so.query <- readRDS("data_integrated_harmony_DUBStepR_2.rds")
cluster.UMAP(so = so.query, group.by = "seurat_clusters") + theme_void() + labs(title = "Elkrewi 2023",subtitle = "Artemia Female Reproductive System")
mc.list <- multiCluster(object = so.query, resolutions = c(0.05, 0.1, 0.2, 0.3, 0.4, 0.5, 0.75, 1, 1.25, 1.5, 1.75, 2), assay = NULL, nworkers = 4, pca_var = 0.9, group_singletons = F, algorithm = 1,return_object = F)
plt.umap_by_cluster <- mc.list$plots
so.query <- mc.list$object
cr_names <- mc.list$resolution_names
cluster.name <- mc.list$cluster_names
assay.pattern <- mc.list$assay_pattern
# visualize cluster configurations
cowplot::plot_grid(plotlist = lapply(plt.umap_by_cluster, function(x) {
    x + theme_void() + labs(title = NULL) + theme(legend.position = "none", plot.subtitle = element_text(hjust = 0.5))
}), ncol = 4)
ms.list <- multiSpecificity(object = so.query, cluster_names = cluster.name, features = NULL, deg_prefilter = T,
    cdi_bins = seq(0, 1, by = 0.01), min.pct = 0.1, n.workers = 16, return_dotplot = T, verbose = T)
df.summary <- ms.list$specificity_summary
df.raw <- ms.list$specificity_raw
# plt.specificity.umap <- ms.list$umap_plot
plt.clust.spec <- ms.list$auc_plot
plt.auc.spec <- ms.list$resolution_plot
plt.auc.dot <- ms.list$dot_plot
max.auc = max(df.summary$auc)
speak <- df.summary$res[which.max(df.summary$auc)]
cowplot::plot_grid(plt.auc.spec + geom_hline(yintercept = max.auc, linetype = "dashed", color = "tomato") +
    geom_vline(xintercept = as.numeric(speak), linetype = "dashed", color = "tomato"), plt.clust.spec,
    plt.auc.dot$`0.5` + theme(legend.position = "right", axis.text.x = element_text(angle = 45,
        hjust = 1)), nrow = 1, rel_widths = c(1, 1.25, 1.25), labels = "AUTO")
cowplot::plot_grid(plotlist = lapply(plt.umap_by_cluster, function(x) {
    x + theme_void() + labs(title = NULL) + theme(legend.position = "none", plot.subtitle = element_text(hjust = 0.5))
}), ncol = 4)
png("umap_by_cluster.png",height=500,width=1200)
cowplot::plot_grid(plotlist = lapply(plt.umap_by_cluster, function(x) {
    x + theme_void() + labs(title = NULL) + theme(legend.position = "none", plot.subtitle = element_text(hjust = 0.5))
}), ncol = 4)
dev.off()
png("cluster_specificity.png",height=500,width=1200)
cowplot::plot_grid(plt.auc.spec + geom_hline(yintercept = max.auc, linetype = "dashed", color = "tomato") +
    geom_vline(xintercept = as.numeric(speak), linetype = "dashed", color = "tomato"), plt.clust.spec,
    plt.auc.dot$`0.2` + theme(legend.position = "right", axis.text.x = element_text(angle = 45,
        hjust = 1)), nrow = 1, rel_widths = c(1, 1.25, 1.25), labels = "AUTO")
dev.off()
png("cluster_specificity_0.1.png",height=500,width=1200)
cowplot::plot_grid(plt.auc.spec + geom_hline(yintercept = max.auc, linetype = "dashed", color = "tomato") +
    geom_vline(xintercept = as.numeric(speak), linetype = "dashed", color = "tomato"), plt.clust.spec,
    plt.auc.dot$`0.1` + theme(legend.position = "right", axis.text.x = element_text(angle = 45,
        hjust = 1)), nrow = 1, rel_widths = c(1, 1.25, 1.25), labels = "AUTO")
dev.off()
#savehistory(file = ".Rhistory_specificity_analysis_plot")
