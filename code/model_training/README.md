# miRBind CNN model retraining

Train the Hejret 2023 CNN model on a subsets of [Manakov 1:1 train dataset](https://zenodo.org/records/14501607) with different sizes.

The subset sizes are:
- 200
- 2 000
- 7 720 = 3 860 (Hejret train set positives) * 2
- 20 000
- 200 000

To run the whole analysis presented in the paper, run 

`bash model_training_analysis.sh -i <Manakov dataset>`

with the Manakov 1:1 train dataset as input.

To retrain on whole Manakov train set or whole Hejret train set, run only steps 2) and 3) from the described pipeline.

## Pipeline description

### 1) Dataset subsampling

Run [make subset](make_subset.py) to create a subset of the dataset containing N positive and N negative samples.

`python make_subset.py --N <n> --dataset <dataset> --output <output>`

### 2) Data encoding

Run the [encode dataset](encode_dataset.py) script to transform dataset in tsv into numpy miRNA x target binding matrix.

`python encode_dataset.py -i <input_dataset> -o <output_prefix>`

### 3) Train the model

Have a look at the [training](training.py) script with model training. The script loads an encoded training dataset and trains the [Hejret 2023 CNN model](https://www.nature.com/articles/s41598-023-49757-z).

Input parameters are:
- *data*: file with the encoded dataset, output from the step (2)
- *labels*: file with the encoded labels, output from (2)
- *dataset_size*: number of samples in the encoded dataset
- *ratio*: ratio of positives to negatives in the train dataset, 10 ~ you have 10x more negatives than positives
- *model*: filename for saving the trained model

`python training.py --data <endoced_dataset> --labels <encoded_labels> --dataset_size <number_of_dataset_samples> --ratio <neg:pos ratio> --model <model_name>`

The trained models are available at https://drive.google.com/drive/folders/1jJQY8XcSYGhS_P7L0MG02u9wM0IGHhle?usp=sharing
