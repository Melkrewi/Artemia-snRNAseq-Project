# Removing background using cellbender
```
module load anaconda3/2022.05
source /mnt/nfs/clustersw/Debian/bullseye/anaconda3/2022.05/activate_anaconda3_2022.05.txt
conda activate CellBender

cellbender remove-background \
     --input  /nfs/scistore18/vicosgrp/melkrewi/Project_snRNA_ovaries/6.FR_gtf_full_data/Afran_1/outs/raw_feature_bc_matrix/ \
     --output ./output \
     --expected-cells 15071 

conda deactivate
```
The file output_filtered.h5 was used for the rest of the analysis
