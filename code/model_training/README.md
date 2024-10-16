# miRBind CNN model retraining

Train miRBind CNN model on a subsets of [Manakov 1:1 train dataset](https://zenodo.org/records/13909173) with different sizes

The subset sizes are:
- 200
- 2 000
- 7 720 = 3 860 (Hejret train set positives) * 2
- 20 000
- 200 000
- 2 524 246 = full Manakov train

## 1) Dataset subsampling

Run [make subset](make_subset.py) to create a subset of the dataset containing N positive and N negative samples.

`python make_subset.py --N <n> --dataset <dataset> --output <output>`

## 2) Data encoding

Run the [encode dataset](encode_dataset.py) script to transform dataset in tsv into numpy miRNA x target binding matrix.

`python encode_dataset.py -i <input_dataset> -o <output_prefix>`

## 3) Train the model

Have a look at the [training](training.py) script with model training.

`python training.py --data <endoced_dataset> --labels <encoded_labels> --dataset_size <number_of_dataset_samples> --ratio <neg:pos ratio> --model <model_name>`

The trained models are available at https://drive.google.com/drive/folders/1jJQY8XcSYGhS_P7L0MG02u9wM0IGHhle?usp=sharing
