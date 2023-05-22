# installation conda


### scvi-tools
```
module load anaconda3/2023.04
source /mnt/nfs/clustersw/Debian/bullseye/anaconda3/2023.04/activate_anaconda3_2023.04.txt
conda create -n scvi-env python=3.9
conda activate scvi-env
conda install pytorch torchvision torchaudio pytorch-cuda=11.7 -c pytorch -c nvidia
conda install jax jaxlib -c conda-forge
conda install scvi-tools -c conda-forge
```
### cellbender
```
module load anaconda3/2023.04

source /mnt/nfs/clustersw/Debian/bullseye/anaconda3/2023.04/activate_anaconda3_2023.04.txt

conda create -n cellBender_reka_test python=3.7

conda activate cellBender_reka_test

conda install -c anaconda pytables

conda install pytorch torchvision torchaudio pytorch-cuda=11.7 -c pytorch -c nvidia

git clone https://github.com/broadinstitute/CellBender.git

pip install -e CellBender
```
### scAR
```
module load anaconda3/2023.04
source /mnt/nfs/clustersw/Debian/bullseye/anaconda3/2023.04/activate_anaconda3_2023.04.txt
git clone https://github.com/Novartis/scar.git
cd scar/
conda env create -f scar-gpu.yml
conda activate scar_gpu
```

