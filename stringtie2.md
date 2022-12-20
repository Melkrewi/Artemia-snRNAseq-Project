# Generating a gtf file for cellranger mkref


### Aligning raw reads to the A. fran genome and generating individual gtfs
I used the adult tissue fastq files: 60541,60542,60543,60544,60545,60546,60547,60548
```
module load samtools
module load hisat2
module load stringtie

prefix="XXXXX"
hisat2 --phred33 -p 50 --novel-splicesite-outfile hisat2/${prefix}_splicesite.txt -S hisat2/${prefix}_accepted_hits.sam -x genome_index -1 ${prefix}_1.fastq -2 ${prefix}_2.fastq --rna-strandness RF
samtools view -bS -o hisat2/${prefix}_accepted_hits.bam hisat2/${prefix}_accepted_hits.sam
samtools sort -o hisat2/${prefix}_accepted_hits.sorted.bam hisat2/${prefix}_accepted_hits.bam
stringtie hisat2/${prefix}_accepted_hits.sorted.bam -o ${prefix}_transcripts.gtf
```
### merging gtf files and keeping transcripts with TPM above 1 and removing the couple unstranded ones:
```
module load stringtie
stringtie --merge 60541_transcripts.gtf 60542_transcripts.gtf 60543_transcripts.gtf 60544_transcripts.gtf 60545_transcripts.gtf 60546_transcripts.gtf 60547_transcripts.gtf 60548_transcripts.gtf -o merged_fran_adult.gtf -F 0 -T 1 #--rf
awk '$7 == "." { next } { print }' merged_fran_adult.gtf > merged_fran_adult_2.gtf
```
### to get a fasta file of transcripts
```
module load cufflinks
gffread -w merged_fran_adult_2.fa -g CHRR_integrated_fran.fa merged_fran_adult_2.gtf
```
