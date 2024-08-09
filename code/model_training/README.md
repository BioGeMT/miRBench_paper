 # miRBind model training
 
## 1) Data encoding

Run the `encode_dataset.py` script to transform dataset in tsv into numpy miRNA x target binding matrix.

The tsv file is expected to have `gene`, `noncodingRNA` and `label` columns.

## 2) Train the model

Have a look at the `training.ipynb` notebook with model training.

### Notes: Dataset sizes and times to encode:

**1:1**

Number of training samples:  2 524 246

Encoding time: 928 s ~ 15 min

Training cca 1,5 hours

**1:10**

Number of training samples:  13 883 353

Encoding time: 5019 s ~ 84 min

Training cca 10 hours

**1:100**

Number of training samples: 127 476 014 !! the set bigger than should be?

Encoding time: 48725 s ~ 812 min ~ 13,5 h