#!/bin/bash
#SBATCH --account=ssamm10
#SBATCH --job-name=Run_enaBrowserTools
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=1
#SBATCH --output=Run_enaBrowserTools.out
#SBATCH --error=Run_enaBrowserTools.err
#SBATCH --exclude=targetid-[05]

# Assign the first argument to the variable download_destination_dir and the second argument to the variable geo_id 
download_destination_dir="$1"
geo_id="$2"

# Get list of SRXs for given geo_id using the GEOparse tool via the python script getSRX_geoparse.py
srx_list=$(python getSRX_geoparse.py "$geo_id")

# ONCE BIOGEMT IS MOUNTED ON THE COMPUTATIONAL SUBNODES!!!, TRY TO RUN LIKE THIS FIRST, IF IT DOESN'T WORK, ADD INSTRUCTIONS TO MANUALLY ADD ALIAS & DIR TO PATH ON .bashrc
current_dir=$(dirname "$0")
enaDataGet_dir="$current_dir/enaBrowserTools/python3/enaDataGet"

# Loop over list of SRXs and download the data using the enaDataGet command
for srx in $srx_list
do
    echo "Downloading $srx" 
    $enaDataGet_dir -f fastq -d "$download_destination_dir" "$srx" 
done
