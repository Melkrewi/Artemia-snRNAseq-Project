# Generating a gtf file for cellranger mkref


### Aligning raw reads to the A. fran genome and generating individual gtfs
I used the adult tissue fastq files: 60541,60542,60543,60544,60545,60546,60547,60548
```
module load cellranger


cellranger mkref --genome=refgenome4cellranger2 --fasta=CHRR_integrated_fran.fa --genes=merged_fran_adult_2.gtf

FASTQS="/nfs/scistore18/vicosgrp/melkrewi/Project_snRNA_ovaries/H32FGDSX5_2_R14367_20221027/demultiplexed/205916/"

cellranger count --transcriptome=/nfs/scistore18/vicosgrp/melkrewi/Project_snRNA_ovaries/6.FR_gtf_full_data/refgenome4cellranger2/ --id=Afran_1 --fastqs=$FASTQS --sample=205916 --include-introns #--localcores=60 --localmem=300 --include-introns

#cellranger reanalyze --id=my_reanalysis --force-cells=8000 --matrix=/nfs/scistore18/vicosgrp/melkrewi/Project_snRNA_ovaries/6.FR_gtf_full_data/Afran_1/outs/raw_feature_bc_matrix.h5
```
