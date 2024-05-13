#INSTRUCTIONS 

git clone miRBench

cd scripts_to_download_eCLIP_Manakov_2022

git clone git@github.com:enasequence/enaBrowserTools.git into this repo

conda env create -f env.yaml

conda activate Run_enaBrowserTools

bash Run_enaBrowserTools.sh <download_destination_dir> <geo_id>

or using slurm, submit your job to the cluster via:

sbatch Run_enaBrowserTools.sh <download_destination_dir> <geo_id>
e.g. sbatch Run_enaBrowserTools.sh /data/biogemt/miRNA_Binding_Site_Prediction/raw_data/eCLIP_Manakov_2022/downloaded_files GSE198250

-------------------
