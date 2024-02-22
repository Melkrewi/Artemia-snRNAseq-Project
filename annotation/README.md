# Generating a gtf file for cellranger mkref


### Aligning raw reads to the A. fran genome and generating individual gtfs
I aligned the RNAseq fastq files of the adult tissue (males and females heads and gonads) to the genome:
```
module load tophat
module load cufflinks
module load bowtie2/2.4.4
module load python/2.7
#run commands on SLURM's srun
module load samtools
module load hisat2
module load stringtie
export TMPDIR=/nfs/scistore18/vicosgrp/melkrewi/Project_snRNA_ovaries_with_W/1.stringtie/
mkdir hisat2
hisat2-build asm_np_female_mkf02_01_09_2023_renamed_final_fin_wscaff.fa genome_index
for i in *_1.fastq
do
    prefix=$(basename $i _1.fastq)
    hisat2 --phred33 -p 50 --novel-splicesite-outfile hisat2/${prefix}_splicesite.txt -S hisat2/${prefix}_accepted_hits.sam -x genome_index -1 ${prefix}_1.fastq -2 ${prefix}_2.fastq --rna-strandness RF
    samtools view -@ 25 -bS -o hisat2/${prefix}_accepted_hits.bam hisat2/${prefix}_accepted_hits.sam
    samtools sort -@ 25 -o hisat2/${prefix}_accepted_hits.sorted.bam hisat2/${prefix}_accepted_hits.bam
    stringtie hisat2/${prefix}_accepted_hits.sorted.bam -o ${prefix}_transcripts.gtf
done
```
### merging gtf files and keeping transcripts with TPM above 1 and removing the couple unstranded ones:
```
module load stringtie
stringtie --merge 60541_transcripts.gtf 60542_transcripts.gtf 60543_transcripts.gtf 60544_transcripts.gtf 60545_transcripts.gtf 60546_transcripts.gtf 60547_transcripts.gtf 60548_transcripts.gtf -o afran_genome_annotation_0.5.gtf -F 0 -T 0.5
awk '$7 == "." { next } { print }' afran_genome_annotation_0.5.gtf > afran_genome_annotation_0.5_2.gtf
```
### to get a fasta file of transcripts
```
module load cufflinks
gffread -w afran_genome_annotation_transcripts_0.5.fa -g asm_np_female_mkf02_01_09_2023_renamed_final_fin_wscaff.fa afran_genome_annotation_0.5_2.gtf
```
