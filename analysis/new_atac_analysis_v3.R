library(harmony)
library(Seurat)
library(ggpubr)
library(dplyr)
library(ggplot2)
library(cowplot)
library(RColorBrewer)
library(corrplot)
library(data.table)
library(Signac)
library(scCustomize)
set.seed(1234)
library(GenomicRanges)
peaks_replicate_3 <- read.table(
  file = "peaks_replicate_3_clean.bed",
  col.names = c("chr", "start", "end")
)
peaks_replicate_4 <- read.table(
  file = "peaks_replicate_4_clean.bed",
  col.names = c("chr", "start", "end")
)
peaks.3 <- makeGRangesFromDataFrame(peaks_replicate_3)
peaks.4 <- makeGRangesFromDataFrame(peaks_replicate_4)
combined.peaks <- reduce(x = c(peaks.3,peaks.4))
peakwidths <- width(combined.peaks)
combined.peaks <- combined.peaks[peakwidths  < 10000 & peakwidths > 20]
data <- readRDS("../data_integrated_harmony_DUBStepR_3_2.rds")
pdf("atac_plots_harmony_10.pdf",width=14,height=6)
subset_3 <- subset(x = data, subset = replicates == "replicate_3")

frag.file <- "~/4.replicate_3/Afran_ATAC/outs/atac_fragments.tsv.gz"

fragcounts <- CountFragments(fragments = frag.file)

atac.frags <- CreateFragmentObject(path = frag.file, cells = colnames(subset_3@assays$RNA))

gtf <- rtracklayer::import('~/4.replicate_3/afran_genome/genes/genes.gtf.gz')
gtf$gene_biotype <- 'protein_coding'
gtf$gene_name <- gtf$gene_id
gtf$tx_id <- gtf$transcript_id
counts <- FeatureMatrix(
  fragments = atac.frags,
  features = combined.peaks,
  cells = colnames(subset_3@assays$RNA)
)
#write.csv(counts,"atac_counts_multiome_1_rpca_4.csv")
atac.assay <- CreateChromatinAssay(
  counts = counts,
  fragments = atac.frags,min.features = 1, annotation=gtf
)

subset_3 <- subset_3[, which((colnames(subset_3) %in% colnames(atac.assay))==TRUE)]
subset_3[["peaks"]] <- atac.assay
DefaultAssay(subset_3) <- "peaks"
subset_3 <- NucleosomeSignal(subset_3)
subset_3 <- TSSEnrichment(subset_3)
subset_3 <- FindTopFeatures(subset_3, min.cutoff = 'q3')
subset_3 <- RunTFIDF(subset_3,method=3)
subset_3 <- RunSVD(subset_3,tol = 1e-12,n=20)

subset_4 <- subset(x = data, subset = replicates == "replicate_4")

frag.file <- "~/5.replicate_4/Afran_ATAC2/outs/atac_fragments.tsv.gz"

fragcounts <- CountFragments(fragments = frag.file)

atac.frags <- CreateFragmentObject(path = frag.file, cells = colnames(subset_4@assays$RNA))
gtf <- rtracklayer::import('~/5.replicate_4/afran_genome/genes/genes.gtf.gz')
gtf$gene_biotype <- 'protein_coding'
gtf$gene_name <- gtf$gene_id
gtf$tx_id <- gtf$transcript_id
counts <- FeatureMatrix(
  fragments = atac.frags,
  features = combined.peaks,
  cells = colnames(subset_4@assays$RNA)
)
#write.csv(counts,"atac_counts_multiome_1_rpca_4.csv")
atac.assay <- CreateChromatinAssay(
  counts = counts,
  fragments = atac.frags,min.features = 1, annotation=gtf
)

subset_4 <- subset_4[, which((colnames(subset_4) %in% colnames(atac.assay))==TRUE)]
subset_4[["peaks"]] <- atac.assay
DefaultAssay(subset_4) <- "peaks"
subset_4 <- NucleosomeSignal(subset_4)
subset_4 <- TSSEnrichment(subset_4)
subset_4 <- FindTopFeatures(subset_4, min.cutoff = 'q3')
subset_4 <- RunTFIDF(subset_4,method=3)
subset_4 <- RunSVD(subset_4,tol = 1e-12,n=20)
combined <- merge(subset_3,subset_4)
combined <- FindTopFeatures(combined, min.cutoff = 10)
combined <- RunTFIDF(combined,method=3)
combined <- RunSVD(combined,tol = 1e-12,n=20)
combined <- RunUMAP(combined, reduction = "lsi", dims = 2:20)
integration.anchors <- FindIntegrationAnchors(
  object.list = list(subset_3,subset_4),
  anchor.features = rownames(subset_4),
  reduction = "rlsi",
  dims = 2:20
)
integrated <- IntegrateEmbeddings(
  anchorset = integration.anchors,
  reductions = combined[["lsi"]],
  new.reduction.name = "integrated_lsi",
  dims.to.integrate = 1:20
)
integrated <- RunUMAP(integrated, reduction = "integrated_lsi", dims = 2:20,min.dist = 0,spread = 4)
subset <- subset(x = data, subset = (replicates == "replicate_4")|(replicates == "replicate_3"))
subset <- subset[, which((colnames(subset) %in% colnames(integrated))==TRUE)]
integrated[['umap.rna']] <- subset[['umap']]
integrated[['umap.atac']] <- integrated[['umap']]
integrated[['umap']] <- integrated[['umap.rna']]
saveRDS(integrated, "ATAC_integrated.rds")
#subset <- subset(x = data, subset = (replicates == "replicate_4")|(replicates == "replicate_3"))
#subset <- subset[, which((colnames(subset) %in% colnames(integrated))==TRUE)]
#subset[['ATAC']] <- integrated[['ATAC']]
#subset[['peaks']] <- integrated[['peaks']]
#subset[['umap.atac']] <- integrated[['umap']]
p1 <- DimPlot(integrated, reduction = "umap", label = TRUE, label.size = 2.5, repel = TRUE) + ggtitle("RNA")
p2 <- DimPlot(integrated, reduction = "umap.atac", label = TRUE, label.size = 2.5, repel = FALSE) + ggtitle("ATAC")
p1 +p2
pdf("atac_umap_final.pdf",width=10,height=5)
p1+p2
VlnPlot(integrated,features="nCount_peaks",split.by="replicates")
dev.off()
library(scrattch.io)
write_dgCMatrix_csv(integrated@assays$peaks@counts, "ATAC_peaks_counts_DUBStep_v3.csv")
#saveRDS(subset, "data_integrated_ATAC_DUBStepR_15_03_2024_v2.rds")
