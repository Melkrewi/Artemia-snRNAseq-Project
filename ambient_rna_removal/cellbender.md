# Removing background using cellbender
### replicate_1
```
module load anaconda3/2022.05
source ~/anaconda3/2022.05/activate_anaconda3_2022.05.txt
conda activate CellBender_cuda
cellbender remove-background \
     --input ~/2.replicate_1/Afran_1/outs/raw_feature_bc_matrix/ \
     --output ./output \
     --expected-cells 18730 --fpr 0.01 0.05 0.1 0.3 0.4 0.5 --cuda --low-count-threshold 2500 --empty-drop-training-fraction 0.4
conda deactivate

```
The file output_filtered.h5 was used for the rest of the analysis
### replicate_2
```
module load anaconda3/2022.05
source ~/anaconda3/2022.05/activate_anaconda3_2022.05.txt
conda activate CellBender_cuda
cellbender remove-background \
     --input  ~/3.replicate_2/Afran_2/outs/raw_feature_bc_matrix/ \
     --output ./output \
     --expected-cells 10000 \
     --total-droplets-included 15000 \
     --epochs 150 --cuda --fpr 0.01 0.05 0.1 0.3 0.35 0.4 0.5 --low-count-threshold 1500 --empty-drop-training-fraction 0.3
conda deactivate
```
### replicate_3
```
module load anaconda3/2022.05
source ~/anaconda3/2022.05/activate_anaconda3_2022.05.txt
conda activate CellBender_cuda
cellbender remove-background \
     --input raw_feature_bc_matrix_rna_only.h5ad \
     --output ./output \
     --expected-cells 10000 \
     --epochs 300 --cuda --low-count-threshold 1000 --empty-drop-training-fraction 0.3 --fpr 0.01 0.05 0.1 0.3 0.35 0.4 0.5
#--learning-rate 0.00005 --z-dim 50 --z-layers 200 --empty-drop-training-fraction 0.5 #--low-count-threshold 5 
conda deactivate
```
### replicate_4
```
module load anaconda3/2022.05
source ~/anaconda3/2022.05/activate_anaconda3_2022.05.txt
conda activate CellBender_cuda
cellbender remove-background \
     --input raw_feature_bc_matrix_rna_only.h5ad \
     --output ./output \
     --expected-cells 10000 \
     --epochs 150 --cuda --low-count-threshold 600 --empty-drop-training-fraction 0.3 --fpr 0.2
conda deactivate
```
