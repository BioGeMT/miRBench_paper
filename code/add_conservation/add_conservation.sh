#!/bin/bash

#SBATCH --account=ssamm10
#SBATCH --job-name=HD_postprocess_on_concat
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=8

python add_conv_scores.py --ifile /data/biogemt/miRNA_Binding_Site_Prediction/miRBench_paper/code/add_conservation/example_input/data.tsv --ofile /data/biogemt/miRNA_Binding_Site_Prediction/miRBench_paper/code/add_conservation/example_output/data_with_conservation.tsv --phyloP_path /data/biogemt/miRNA_Binding_Site_Prediction/conversation_bw/hg38.phyloP100way.bw --phastCons_path /data/biogemt/miRNA_Binding_Site_Prediction/conversation_bw/hg38.phastCons100way.bw