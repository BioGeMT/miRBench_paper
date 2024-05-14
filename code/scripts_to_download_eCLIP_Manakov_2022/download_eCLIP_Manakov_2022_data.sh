#!/bin/bash
#SBATCH --account=ssamm10
#SBATCH --job-name=download_eCLIP_Manakov_2022_data
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=1
#SBATCH --output=download_eCLIP_Manakov_2022_data.out
#SBATCH --error=download_eCLIP_Manakov_2022_data.err
#SBATCH --exclude=targetid-[05]

# Assign the first argument to the variable download_destination_dir and the second argument to the variable geo_id 
dest_dir="$1"
geo_id="$2"

# Get list of SRXs for given geo_id using the GEOparse tool via the python script getSRX_geoparse.py
srx_list=$(python getSRX_geoparse.py "$geo_id" "$dest_dir")

# Get the current directory
if [ -n "$SLURM_JOB_ID" ]; then
    current_dir=${SLURM_SUBMIT_DIR}
else
    current_dir=$(dirname "$0")
fi 

# Assign the path to the enaDataGet command via the current directory
enaDataGet_dir="$current_dir/enaBrowserTools/python3/enaDataGet"

# Loop over list of SRXs and download the data using the enaDataGet command
for srx in $srx_list
do
    echo "Downloading $srx" 
    $enaDataGet_dir -f fastq -d "$dest_dir" "$srx" 
done

# Print message to indicate date and time when all data was downloaded
echo "All eCLIP Manakov 2022 data was downloaded on " `date` " by Stephanie Sammut."