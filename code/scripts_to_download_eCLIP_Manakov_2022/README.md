# Description

This package downloads all the runs (SRR) of all the experiments (SRX) of all the samples (GSM) under a given GEO series / GEO ID (GSE) number. 

# Commands to run on command line

`git clone git@github.com:enasequence/enaBrowserTools.git`

`conda env create -f env.yaml`

`conda activate env_to_download_eCLIP_Manakov_2022_data`

`bash download_eCLIP_Manakov_2022_data.sh <download_destination_dir> <geo_id>`

or using slurm, edit SBATCH directives in .sh and submit your job to the cluster via:

`sbatch download_eCLIP_Manakov_2022_data.sh <download_destination_dir> <geo_id>`

e.g. `sbatch download_eCLIP_Manakov_2022_data.sh /data/biogemt/miRNA_Binding_Site_Prediction/raw_data/eCLIP_Manakov_2022 GSE198250`
