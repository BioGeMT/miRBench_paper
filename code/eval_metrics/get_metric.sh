#!/bin/bash

#SBATCH --account=ssamm10
#SBATCH --job-name=get_metric
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=8

## declare an array variable
declare -a datasets=("eCLIP_Klimentova2022_test") # "CLASH_Hejret2023_test" "eCLIP_Manakov2022_test" "eCLIP_Manakov2022_leftout")

## declare an array variable
declare -a metrics=("auc_pr" "auc_roc" "avg_p_score")

# now loop through the above array
for dataset in "${datasets[@]}"; do
    for metric in "${metrics[@]}"; do
        python ./code/eval_metrics/get_metric.py \
            --ifile ./code/eval_metrics/example_input/AGO2_${dataset}_predictions.tsv \
            --ofile ./code/eval_metrics/example_output/AGO2_${dataset}_${metric}.tsv \
            --metric ${metric}
    done
done
