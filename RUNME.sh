#!/bin/bash

## This file can be executed and will reproduce the methods of the paper.

# Step 1: Data Preprocessing

## Installation of required software etc

### TODO: Add the code to install the required software

## Download the data (Chimeric RNA dataset)

### All information for downloading the data is in:
### /code/scripts_to_download_eCLIP_Manakov_2022/
### We will not reproduce the download here, but we will assume that the data is in the data/ directory

## Process the data

### TODO: Add the code to process the data using HybriDetector
### probably here starts a loop over the datasets and everything following is in it
### keep the list of dataset explicit in each loop

# dummy_code
# for dataset in samplex, sampley, samplez, ...
# do
#     download_data $dataset > data/$dataset/raw_data.txt
# done
# dummy_code

# From HybriDetector output to list of miRNA-target interactions

### TODO: Add the code to convert the output of HybriDetector to a list of miRNA-target interactions
### miRNAfinder (how was it run for each dataset)
### keep each argument in a separate line for better readability

# python miRNAfinder/miRNAfinder.py \
#     --ifile data/samplex/HybriDetector_output.tsv \
#     --ofile data/samplex/miRNA_target_interactions.tab

## Produce negative samples
### TODO: Add the code to produce negative samples

# dummy_code
# python dev/negative_samples.py \
#     --ifile data/samplex/miRNA_target_interactions.tab \
#     --ofile data/samplex/miRNA_target_interactions.posneg.tab \
#     --negratio 100

## At this point, we have the positive and negative samples for each dataset in a separate file for each dataset
## for ease of use we will move all the files to a single directory, and change their name to include the sample name

# dummy_code
# for dataset in samplex, sampley, samplez, ...
# do
#     mv \
#     data/$dataset/miRNA_target_interactions.posneg.tab \
#     data/posneg_samples/$dataset.miRNA_target_interactions.posneg.tab
# done
# dummy_code

## Train/Test Split (this uses the specific column produced by miRNAfinder OR a random split)

# dummy_code
# for posneg_dataset in data/posneg_samples/*.miRNA_target_interactions.posneg.tab
# do
#     python dev/train_test_split.py \
#         --ifile $posneg_dataset \
#         --trainfile data/train_test_splits/$(basename $posneg_dataset .tab).train.tab \
#         --testfile data/train_test_splits/$(basename $posneg_dataset .tab).test.tab \
#         --trainratio hardcoded # or -trainratio 0.8
# done
# dummy_code

## At this point, we have the train and test splits for each dataset in a separate file for each dataset
## These files can now be used to train or evaluate models

# Step 2: Model Training

...

# Step 3: Model Evaluation

## Here goes code for each model that will produce a file foreach 
## data/train_test_splits/*.test.tab with a column called "pred_score" with descending scores
## i.e. HIGH SCORE IS ALWAYS BETTER

## Evaluation of each model on each dataset

# dummy_code
# for dataset in data/train_test_splits/*.test.tab
# do
#     for model in model1, model2, model3, ...
#     do
#         python dev/evaluate_model.py \
#             --ifile $dataset \
#             --ofile data/evaluation_results/$(basename $dataset .tab).eval.tab \
#             --model $model_name \
#             --model_predictions results/model_predictions/$model_name/$(basename $dataset .tab).pred.tab
#     done
# done

# We now join all the prediction results per dataset as columns in a final table
# The rows are in the same order, the column name is the model name

# dummy_code
# for dataset in samplex, sampley, samplez, ...
# do
#     python dev/join_model_predictions.py \
#         --ifolder model_predictions \
#         --dataset $dataset \
#         --ofile results/prediction/$dataset.eval.all_models.tab
# done
# dummy_code

# Step 4: Model Comparison

## Precision-Recall curves

# dummy_code
# for dataset in samplex, sampley, samplez, ...
# do
#     python dev/plot_precision_recall.py \
#         --ifile results/prediction/$dataset.eval.all_models.tab \
#         --oplot results/plots/$dataset.precision_recall.png \
#         --ofile results/evaluation/$dataset.precision_recall.tab
# done
# dummy_code

## ROC curves

# dummy_code
# for dataset in samplex, sampley, samplez, ...
# do
#     python dev/plot_roc.py \
#         --ifile results/prediction/$dataset.eval.all_models.tab \
#         --oplot results/plots/$dataset.roc.png \
#         --ofile results/evaluation/$dataset.roc.tab
# done
# dummy_code

# ## Boxplots of PR-AUC

# dummy_code
# python dev/PR_AUC_boxplot.py \
#     --ifolder results/evaluation \
#     --oplot results/plots/boxplot_pr_auc.png \
#     --ofile results/evaluation/boxplot_pr_auc.tab
# dummy_code


