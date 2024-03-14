
# Running velocyto for Drosophila
### Downloading the reads
```
wget ftp://ftp.ebi.ac.uk/pub/databases/microarray/data/experiment/MTAB/E-MTAB-10519/FCA62_Female_ovary_adult_5dWT_Nystul_All_Nuclei_S57_L004_R1_001.fastq.gz
wget ftp://ftp.ebi.ac.uk/pub/databases/microarray/data/experiment/MTAB/E-MTAB-10519/FCA62_Female_ovary_adult_5dWT_Nystul_All_Nuclei_S57_L004_R2_001.fastq.gz
wget ftp://ftp.ebi.ac.uk/pub/databases/microarray/data/experiment/MTAB/E-MTAB-10519/FCA62_Female_ovary_adult_5dWT_Nystul_All_Nuclei_S58_L004_R1_001.fastq.gz
wget ftp://ftp.ebi.ac.uk/pub/databases/microarray/data/experiment/MTAB/E-MTAB-10519/FCA62_Female_ovary_adult_5dWT_Nystul_All_Nuclei_S58_L004_R2_001.fastq.gz
wget ftp://ftp.ebi.ac.uk/pub/databases/microarray/data/experiment/MTAB/E-MTAB-10519/FCA62_Female_ovary_adult_5dWT_Nystul_All_Nuclei_S62_L004_R1_001.fastq.gz
wget ftp://ftp.ebi.ac.uk/pub/databases/microarray/data/experiment/MTAB/E-MTAB-10519/FCA62_Female_ovary_adult_5dWT_Nystul_All_Nuclei_S62_L004_R2_001.fastq.gz
wget ftp://ftp.ebi.ac.uk/pub/databases/microarray/data/experiment/MTAB/E-MTAB-10519/FCA62_Female_ovary_adult_5dWT_Nystul_All_Nuclei_S63_L004_R1_001.fastq.gz 
wget ftp://ftp.ebi.ac.uk/pub/databases/microarray/data/experiment/MTAB/E-MTAB-10519/FCA62_Female_ovary_adult_5dWT_Nystul_All_Nuclei_S63_L004_R2_001.fastq.gz
wget ftp://ftp.ebi.ac.uk/pub/databases/microarray/data/experiment/MTAB/E-MTAB-10519/FCA63_Female_ovary_adult_5dWT_Nystul_Big_Nuclei_S64_L004_R1_001.fastq.gz
wget ftp://ftp.ebi.ac.uk/pub/databases/microarray/data/experiment/MTAB/E-MTAB-10519/FCA63_Female_ovary_adult_5dWT_Nystul_Big_Nuclei_S64_L004_R2_001.fastq.gz
wget ftp://ftp.ebi.ac.uk/pub/databases/microarray/data/experiment/MTAB/E-MTAB-10519/FCA63_Female_ovary_adult_5dWT_Nystul_Big_Nuclei_S66_L004_R1_001.fastq.gz
wget ftp://ftp.ebi.ac.uk/pub/databases/microarray/data/experiment/MTAB/E-MTAB-10519/FCA63_Female_ovary_adult_5dWT_Nystul_Big_Nuclei_S66_L004_R2_001.fastq.gz
wget ftp://ftp.ebi.ac.uk/pub/databases/microarray/data/experiment/MTAB/E-MTAB-10519/FCA63_Female_ovary_adult_5dWT_Nystul_Big_Nuclei_S71_L004_R1_001.fastq.gz
wget ftp://ftp.ebi.ac.uk/pub/databases/microarray/data/experiment/MTAB/E-MTAB-10519/FCA63_Female_ovary_adult_5dWT_Nystul_Big_Nuclei_S71_L004_R2_001.fastq.gz
wget ftp://ftp.ebi.ac.uk/pub/databases/microarray/data/experiment/MTAB/E-MTAB-10519/FCA63_Female_ovary_adult_5dWT_Nystul_Big_Nuclei_S72_L004_R1_001.fastq.gz
wget ftp://ftp.ebi.ac.uk/pub/databases/microarray/data/experiment/MTAB/E-MTAB-10519/FCA63_Female_ovary_adult_5dWT_Nystul_Big_Nuclei_S72_L004_R2_001.fastq.gz
```
### cellranger
```
module load cellranger


cellranger mkref --genome=refgenome4cellranger2 --fasta=Drosophila_melanogaster.BDGP6.32.dna.toplevel.fa --genes=Drosophila_melanogaster.BDGP6.32.108.chr.gtf

FASTQS="/nfs/scistore18/vicosgrp/melkrewi/Project_snRNA_ovaries/32.reanalyze_drosophila_data/FASTQ/"

cellranger count --transcriptome=refgenome4cellranger2 --id=dmel_1 --fastqs=$FASTQS --sample=FCA62_Female_ovary_adult_5dWT_Nystul_All_Nuclei,FCA63_Female_ovary_adult_5dWT_Nystul_Big_Nuclei --include-introns
```
### Velocyto
```
module load anaconda3/2022.05

source ~/anaconda3/2022.05/activate_anaconda3_2022.05.txt

#cut -f1-2,4- -d";" genes.gtf > genes_2.gtf

conda activate velocyto

module load samtools

velocyto run10x --samtools-threads 40 /nfs/scistore18/vicosgrp/melkrewi/Project_snRNA_ovaries/32.reanalyze_drosophila_data/analysis/dmel_1/ genes.gtf

conda deactivate

```


