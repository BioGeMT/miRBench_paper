#!/bin/bash

#SBATCH --account=ssamm10
#SBATCH --job-name=get_metric
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=8

## declare an array variable
declare -a datasets=("CLASH_Hejret2023" "eCLIP_Klimentova2022" "eCLIP_Manakov2022")

## declare an array variable
declare -a ratios=("1" "10" "100")

# now loop through the above array
for dataset in "${datasets[@]}"; do
    for ratio in "${ratios[@]}"; do
        python /data/biogemt/miRNA_Binding_Site_Prediction/miRBench_paper/code/eval_metrics/get_metric.py \
            --ifile /data/biogemt/miRNA_Binding_Site_Prediction/mirbench/mirbench_out/AGO2_${dataset}_${ratio}_predictions.tsv \
            --ofile /data/biogemt/miRNA_Binding_Site_Prediction/eval_metrics/avg_p_score/AGO2_${dataset}_${ratio}_avg_p_score.tsv \
            --metric avg_p_score
    done
done
