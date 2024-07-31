
# Running velocyto for each replicate
### replicate_1
genes.gtf.gz is from the cellranger output and the Afran_1 folder as well.
```
module load java
export TMPDIR=~/1.replicate_1/
#run commands on SLURM's srun
module load anaconda3/2022.05

source ~/anaconda3/2022.05/activate_anaconda3_2022.05.txt

gunzip genes.gtf.gz
cut -f1-2,4- -d";" genes.gtf > genes_2.gtf

conda activate velocyto

module load samtools

velocyto run10x --samtools-threads 40 ~/2.replicate_1/Afran_1/ genes_2.gtf

conda deactivate
```
### replicate_2
```
module load java
#run commands on SLURM's srun
module load anaconda3/2022.05

source ~/anaconda3/2022.05/activate_anaconda3_2022.05.txt

gunzip genes.gtf.gz
cut -f1-2,4- -d";" genes.gtf > genes_2.gtf

conda activate velocyto

module load samtools

velocyto run10x --samtools-threads 40 ~/3.replicate_2/Afran_2/ genes_2.gtf

conda deactivate

```
### replicate_3
```
module load anaconda3/2022.05

source ~/anaconda3/2022.05/activate_anaconda3_2022.05.txt

gunzip genes.gtf.gz
cut -f1-2,4- -d";" genes.gtf > genes_2.gtf

conda activate velocyto

module load samtools

velocyto run --samtools-threads 40 ~/4.replicate_3/Afran_ATAC/outs/gex_possorted_bam.bam genes_2.gtf

conda deactivate

```
### replicate_4
```
module load anaconda3/2022.05

source ~/anaconda3/2022.05/activate_anaconda3_2022.05.txt

gunzip genes.gtf.gz
cut -f1-2,4- -d";" genes.gtf > genes_2.gtf

conda activate velocyto

module load samtools

velocyto run --samtools-threads 40 ~/5.replicate_4/Afran_ATAC2/outs/gex_possorted_bam.bam genes_2.gtf

conda deactivate
```
