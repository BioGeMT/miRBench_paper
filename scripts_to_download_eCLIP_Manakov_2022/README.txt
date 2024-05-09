#INSTRUCTIONS 

git clone git@github.com:enasequence/enaBrowserTools.git

conda env create -f env.yaml

conda activate Run_enaBrowserTools

bash Run_enaBrowserTools.sh <download_destination_dir> <geo_id>

-------------------
#MY NOTES

manually add to your .bashrc - might not be needed as i run .sh
alias enaDataGet='/homes/ssamm10/test-enaBrowserTools/enaBrowserTools/python3/enaDataGet'

# cd enaBrowserTools (git cloned from git@github.com:enasequence/enaBrowserTools.git)
# Added directory of the enaGetData command to PATH in .bashrc
# Added alias for the enaGetData command in .bashrc (remember to save, exit and apply changes via `source ~/.bashrc`)
# conda activate test (or whatever env contains installation of GEOparse)
# Created a python script that gets the SRXs from the GEOparse tool called getSRX_geoparse.py in enaBrowserTools

# Run `sbatch Run_enaBrowserTools.sh <download_destination_dir> <geo_id> 
# e.g. sbatch Run_enaBrowserTools.sh /data/biogemt/miRNA_Binding_Site_Prediction/raw_data/eCLIP_Manakov_2022/downloaded_files GSE198250`