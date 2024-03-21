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
data <- readRDS("data_integrated_harmony_DUBStepR_3_2.rds")
subset_3 <- subset(x = data, subset = replicates == "replicate_3")
inputdata.10x <- Read10X_h5("/nfs/scistore18/vicosgrp/melkrewi/Project_snRNA_ovaries_with_W/4.replicate_3/Afran_ATAC/outs/raw_feature_bc_matrix.h5")
frag.file <- "/nfs/scistore18/vicosgrp/melkrewi/Project_snRNA_ovaries_with_W/4.replicate_3/Afran_ATAC/outs/atac_fragments.tsv.gz"
atac_counts <- inputdata.10x$Peaks
atac_counts <- atac_counts[, which((colnames(atac_counts) %in% colnames(subset_3@assays$RNA))==TRUE)]
subset_3 <- subset_3[, which((colnames(subset_3) %in% colnames(atac_counts))==TRUE)]
gtf <- rtracklayer::import('/nfs/scistore18/vicosgrp/melkrewi/Project_snRNA_ovaries_with_W/4.replicate_3/afran_genome/genes/genes.gtf.gz')
gtf$gene_biotype <- 'protein_coding'
gtf$gene_name <- gtf$gene_id
gtf$tx_id <- gtf$transcript_id
chrom_assay <- CreateChromatinAssay(
   counts = atac_counts,
   sep = c(":", "-"),
   fragments = frag.file,
   min.cells = 1, annotation=gtf)
subset_3[["ATAC"]] <- chrom_assay
DefaultAssay(subset_3) <- "ATAC"
subset_3 <- NucleosomeSignal(subset_3)
subset_3 <- TSSEnrichment(subset_3)
saveRDS(subset_3, "replicate_3_TSS_enrichment.rds")
peaks <- CallPeaks(subset_3, mac2.path="/nfs/scistore18/vicosgrp/melkrewi/.conda/envs/macs2/bin/macs2",group.by ="seurat_clusters")
rtracklayer::export.bed(peaks, "peaks_replicate_3.bed")
