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
data <- readRDS("data_integrated_harmony_DUBStepR_2.rds")
pdf("atac_plots_harmony_10.pdf",width=14,height=6)
subset_3 <- subset(x = data, subset = replicates == "replicate_4")
inputdata.10x <- Read10X_h5("/nfs/scistore18/vicosgrp/melkrewi/Project_snRNA_ovaries_with_W/5.replicate_4/Afran_ATAC2/outs/raw_feature_bc_matrix.h5")
frag.file <- "/nfs/scistore18/vicosgrp/melkrewi/Project_snRNA_ovaries_with_W/5.replicate_4/Afran_ATAC2/outs/atac_fragments.tsv.gz"
atac_counts <- inputdata.10x$Peaks
atac_counts <- atac_counts[, which((colnames(atac_counts) %in% colnames(subset_3@assays$RNA))==TRUE)]
#atac_counts <- atac_counts[rownames(atac_counts) %like% "LG", ]
#write.csv(atac_counts,"atac_counts_multiome_2_rpca_4.csv")
subset_3 <- subset_3[, which((colnames(subset_3) %in% colnames(atac_counts))==TRUE)]
gtf <- rtracklayer::import('/nfs/scistore18/vicosgrp/melkrewi/Project_snRNA_ovaries_with_W/5.replicate_4/afran_genome/genes/genes.gtf.gz')
gtf$gene_biotype <- 'protein_coding'
gtf$gene_name <- gtf$gene_id
gtf$tx_id <- gtf$transcript_id
#gtf <- gtf[gtf@seqnames %like% "LG", ]
chrom_assay <- CreateChromatinAssay(
   counts = atac_counts,
   sep = c(":", "-"),
   fragments = frag.file,
   min.cells = 1, annotation=gtf)
subset_3[["ATAC"]] <- chrom_assay
DefaultAssay(subset_3) <- "ATAC"
subset_3 <- NucleosomeSignal(subset_3)
subset_3 <- TSSEnrichment(subset_3)
VlnPlot(
  object = subset_3,
  features = c("nCount_RNA", "nCount_ATAC", "TSS.enrichment","nucleosome_signal"),
  ncol = 4,
  pt.size = 0,group.by='orig.ident'
)
subset_3 <- subset(
  x = subset_3,
  subset = nCount_ATAC < 7000 &
    nCount_RNA < 20000 &
    nCount_ATAC > 100 &
    nCount_RNA > 200 &
    nucleosome_signal < 1.7 & TSS.enrichment > 0.5 &
    TSS.enrichment < 10
)
VlnPlot(
  object = subset_3,
  features = c("nCount_RNA", "nCount_ATAC", "TSS.enrichment","nucleosome_signal"),
  ncol = 4,
  pt.size = 0,group.by='orig.ident'
)
peaks <- CallPeaks(subset_3, mac2.path="/nfs/scistore18/vicosgrp/melkrewi/.conda/envs/macs2/bin/macs2",group.by ="seurat_clusters")
macs2_counts <- FeatureMatrix(
  fragments = Fragments(subset_3),
  features = peaks,
  cells = colnames(subset_3)
)
subset_3[["peaks"]] <- CreateChromatinAssay(
  counts = macs2_counts,
  fragments = frag.file,
  annotation = gtf
)
DefaultAssay(subset_3) <- "peaks"
#subset_3 <- subset(
#  x = subset_3,
#  subset = nCount_ATAC > 50 &
#    nCount_RNA > 1
#)
#subset_3 <- RunTFIDF(subset_3,method=3)
#subset_3 <- FindTopFeatures(subset_3, min.cutoff = 'q5')
#subset_3 <- RunSVD(subset_3,tol = 1e-12)
#DepthCor(subset_3,n=NULL)
#subset_3 <- RunUMAP(subset_3, reduction = 'lsi', dims = 2:50, reduction.name = "umap.atac", reduction.key = "atacUMAP_",min.dist = 0.05,spread = 3)
#p1 <- DimPlot(subset_3, reduction = "umap", label = TRUE, label.size = 2.5, repel = TRUE) + ggtitle("RNA")
#p2 <- DimPlot(subset_3 , reduction = "umap.atac", label = TRUE, label.size = 2.5, repel = FALSE) + ggtitle("ATAC")
#p1 +p2
subset_3 <- RunTFIDF(subset_3,method=3)
subset_3 <- FindTopFeatures(subset_3, min.cutoff = 'q3')
subset_3 <- RunSVD(subset_3,tol = 1e-12,n = 20)
DepthCor(subset_3,n=NULL)
subset_3 <- RunUMAP(subset_3, reduction = 'lsi', dims = c(2,3,4,5,6,7,8,9,10,11,12,13,14,15,16,17,18,19,20), reduction.name = "umap.atac", reduction.key = "atacUMAP_",min.dist = 0,spread = 7)
p1 <- DimPlot(subset_3, reduction = "umap", label = TRUE, label.size = 2.5, repel = TRUE) + ggtitle("RNA")
p2 <- DimPlot(subset_3 , reduction = "umap.atac", label = TRUE, label.size = 2.5, repel = FALSE) + ggtitle("ATAC")
p1 +p2
dev.off()
subset_4 <- subset(x = data, subset = replicates == "replicate_3")

frag.file <- "/nfs/scistore18/vicosgrp/melkrewi/Project_snRNA_ovaries_with_W/4.replicate_3/Afran_ATAC/outs/atac_fragments.tsv.gz"

fragcounts <- CountFragments(fragments = frag.file)

atac.frags <- CreateFragmentObject(path = frag.file, cells = colnames(subset_4@assays$RNA))

counts <- FeatureMatrix(
  fragments = atac.frags,
  features = granges(subset_3),
  cells = colnames(subset_4@assays$RNA)
)
#write.csv(counts,"atac_counts_multiome_1_rpca_4.csv")
atac.assay <- CreateChromatinAssay(
  counts = counts,
  min.features = 20,
  fragments = atac.frags
)

pbmc.atac <- CreateSeuratObject(counts = atac.assay, assay = "peaks")
pbmc.atac <- FindTopFeatures(pbmc.atac, min.cutoff = 'q3')
pbmc.atac <- RunTFIDF(pbmc.atac,method=3)
pbmc.atac <- RunSVD(pbmc.atac,tol = 1e-12,n=20)
pbmc.multi <- subset_3
pbmc.atac$dataset <- "multiome_3"
pbmc.multi$dataset <- "multiome_4"
pbmc.combined <- merge(pbmc.multi,pbmc.atac)
pbmc.combined <- FindTopFeatures(pbmc.combined, min.cutoff = 10)
pbmc.combined <- RunTFIDF(pbmc.combined,method=3)
pbmc.combined <- RunSVD(pbmc.combined,tol = 1e-12,n=20)
pbmc.combined <- RunUMAP(pbmc.combined, reduction = "lsi", dims = 2:20)
p1 <- DimPlot(pbmc.combined, group.by = "dataset")
p1
integration.anchors <- FindIntegrationAnchors(
  object.list = list(pbmc.multi, pbmc.atac),
  anchor.features = rownames(pbmc.multi),
  reduction = "rlsi",
  dims = 2:20
)
integrated <- IntegrateEmbeddings(
  anchorset = integration.anchors,
  reductions = pbmc.combined[["lsi"]],
  new.reduction.name = "integrated_lsi",
  dims.to.integrate = 1:20
)

integrated <- RunUMAP(integrated, reduction = "integrated_lsi", dims = 2:20)
p2
DimPlot(integrated, group.by = "dataset")
integrated <- IntegrateEmbeddings(
  anchorset = integration.anchors,
  reductions = pbmc.combined[["lsi"]],
  new.reduction.name = "integrated_lsi",
  dims.to.integrate = 1:20
)

integrated <- RunUMAP(integrated, reduction = "integrated_lsi", dims = 2:20,min.dist = 0,spread = 7)
subset <- subset(x = data, subset = (replicates == "replicate_4")|(replicates == "replicate_3"))
subset <- subset[, which((colnames(subset) %in% colnames(integrated))==TRUE)]
subset[['ATAC']] <- integrated[['ATAC']]
subset[['umap.atac']] <- integrated[['umap']]
p1 <- DimPlot(subset, reduction = "umap", label = TRUE, label.size = 2.5, repel = TRUE) + ggtitle("RNA")
p2 <- DimPlot(subset , reduction = "umap.atac", label = TRUE, label.size = 2.5, repel = FALSE) + ggtitle("ATAC")
p1 +p2
#subset_clean <- subset[,  subset$seurat_clusters %in% c("12","14","4","9","1","5","6","10","7","15","16","0","3","8")]
#p1 <- DimPlot(subset_clean, reduction = "umap", label = TRUE, label.size = 2.5, repel = TRUE) + ggtitle("RNA")
#p2 <- DimPlot(subset_clean , reduction = "umap.atac", label = TRUE, label.size = 2.5, repel = FALSE) + ggtitle("ATAC")
#p1 +p2
#p1 <- DimPlot(subset_clean, reduction = "umap", label = TRUE, label.size = 2.5, repel = TRUE) + ggtitle("RNA")
#p2 <- DimPlot(subset_clean , reduction = "umap.atac", label = TRUE, label.size = 2.5, repel = FALSE) + ggtitle("ATAC")
#p1 +p2
pdf("atac_umap_final.pdf",width=10,height=5)
p1+p2
dev.off()
saveRDS(subset, "data_integrated_ATAC_DUBStepR.rds")
