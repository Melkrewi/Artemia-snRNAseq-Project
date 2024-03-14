# making a web browser page


### downloading cellbrowser
```
load anaconda3/2022.05
source /mnt/nfs/clustersw/Debian/bullseye/anaconda3/2022.05/activate_anaconda3_2022.05.txt
conda create -n cellbrowser
conda activate cellbrowser
conda install -c bioconda ucsc-cell-browser
```
### generating html page
```
load anaconda3/2022.05
source /mnt/nfs/clustersw/Debian/bullseye/anaconda3/2022.05/activate_anaconda3_2022.05.txt
conda activate cellbrowser
cp /nfs/scistore18/vicosgrp/melkrewi/Project_snRNA_ovaries/22.cellblender/analysis_7/data_filt.h5ad .
module load R/4.2.1
cbImportSeurat -i data_integrated_harmony_DUBStepR_3_2_renamed.rds -o afransinglenucleusatlas
cd afransinglenucleusatlas
cbBuild -o ./public_html/cb
```
