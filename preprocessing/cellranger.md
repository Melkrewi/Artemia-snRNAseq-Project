# Preprocessing data with cellranger

### making reference and getting counts
### replicate_1 (10x RNAseq)
```
module load cellranger/5.0.0

cellranger mkref --genome=refgenome4cellranger2 --fasta=~/asm_np_female_mkf02_01_09_2023_renamed_final_fin_wscaff.fa --genes=~/afran_genome_annotation_0.5_2.gtf

FASTQS="~/H32FGDSX5_2_R14367_20221027/demultiplexed/205916/"

cellranger count --transcriptome=refgenome4cellranger2 --id=Afran_1 --fastqs=$FASTQS --sample=205916 --include-introns
```
### replicate_2 (10x RNAseq)
```
module load cellranger/5.0.0

cellranger mkref --genome=refgenome4cellranger2 --fasta=~/asm_np_female_mkf02_01_09_2023_renamed_final_fin_wscaff.fa --genes=~/afran_genome_annotation_0.5_2.gtf

FASTQS="~/H335JDSX5_1_R14694_20230119/demultiplexed/216143/"

cellranger count --transcriptome=refgenome4cellranger2/ --id=Afran_2 --fastqs=$FASTQS --sample=216143 --include-introns
```
### replicate_3 (10x multiome)
Example of conflig file and libraries.csv is in the folder.
```
~/cellranger-arc-2.0.2/cellranger-arc mkref --config=config_file.txt

~/cellranger-arc-2.0.2/cellranger-arc count --reference=afran_genome --id=Afran_ATAC --libraries=libraries.csv
```
### replicate_4 (10x multiome)
```
module load cellranger/5.0.0

~/cellranger-arc-2.0.2/cellranger-arc mkref --config=config_file.txt

~/cellranger-arc-2.0.2/cellranger-arc count --reference=afran_genome --id=Afran_ATAC2 --libraries=libraries.csv
```


