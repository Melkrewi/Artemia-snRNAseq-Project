# Generating a gtf file for cellranger mkref


## Kmer analysis 
```
prefix="60548"
hisat2 --phred33 -p 50 --novel-splicesite-outfile hisat2/${prefix}_splicesite.txt -S hisat2/${prefix}_accepted_hits.sam -x genome_index -1 ${prefix}_1.fastq -2 ${prefix}_2.fastq --rna-strandness RF
samtools view -bS -o hisat2/${prefix}_accepted_hits.bam hisat2/${prefix}_accepted_hits.sam
samtools sort -o hisat2/${prefix}_accepted_hits.sorted.bam hisat2/${prefix}_accepted_hits.bam
stringtie hisat2/${prefix}_accepted_hits.sorted.bam -o ${prefix}_transcripts.gtf
```
