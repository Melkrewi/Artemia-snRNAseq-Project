### Mapping
We first mapped the expressed transcripts (list of expressed genes included) in our single cell data set to the drosophila CDS (removed isoforms, text file with transcripts included) using the mapping script from SAMAP as shown below. The drosophila cds was downloaded from: https://ftp.flybase.net/genomes/dmel/dmel_r6.31_FB2019_06/fasta/
```
module load blast
/nfs/scistore18/vicosgrp/melkrewi/Project_snRNA_ovaries/8.integrate_with_drosophila/5.integrate/SAMap/map_genes.sh --tr1 afran_genome_annotation_transcripts_0.5_expressed.fa --t1 nucl --n1 AR --tr2 unique_isoforms_dmel-all-CDS-r6.31.fasta --t2 nucl --n2 DM --threads 60
```
### Integration
The code for integration is in the included jupyter notebook.
