#!/usr/bin/env python
# coding: utf-8

import random
import numpy as np
import argparse
import time
import tensorflow as tf
import matplotlib.pyplot as plt
from tensorflow.keras.layers import (
                                BatchNormalization, LeakyReLU,
                                Input, Dense, Conv2D,
                                MaxPooling2D, Flatten, Dropout)
from tensorflow.keras.optimizers import Adam
from tensorflow.keras.utils import Sequence
import os


def make_architecture():
    """
    Build model architecture
    The parameters are based on the CNN model presented in miRBind paper (2022) https://doi.org/10.3390/genes13122323
    Link to the original implementation: https://github.com/ML-Bioinfo-CEITEC/miRBind/blob/main/Additional_scripts/hyperparams.ipynb
    """
    cnn_num = 6
    kernel_size = 5
    pool_size = 2
    dropout_rate = 0.3
    dense_num = 2

    x = Input(shape=(50,20,1),
                       dtype='float32', name='main_input'
                       )
    main_input = x

    for cnn_i in range(cnn_num):
        x = Conv2D(
            filters=32 * (cnn_i + 1),
            kernel_size=(kernel_size, kernel_size),
            padding="same",
            data_format="channels_last",
            name="conv_" + str(cnn_i + 1))(x)
        x = LeakyReLU()(x)
        x = BatchNormalization()(x)
        x = MaxPooling2D(pool_size=(pool_size, pool_size), padding='same', name='Max_' + str(cnn_i + 1))(x)
        x = Dropout(rate=dropout_rate)(x)

    x = Flatten(name='2d_matrix')(x)

    for dense_i in range(dense_num):
        neurons = 32 * (cnn_num - dense_i)
        x = Dense(neurons)(x)
        x = LeakyReLU()(x)
        x = BatchNormalization()(x)
        x = Dropout(rate=dropout_rate)(x)

    main_output = Dense(1, activation='sigmoid', name='main_output')(x)

    model = tf.keras.Model(inputs=[main_input], outputs=[main_output], name='arch_00')
    
    return model


def compile_model():
    tf.keras.backend.clear_session()
    model = make_architecture()
    
    opt = Adam(
        learning_rate=0.00152,
        beta_1=0.9,
        beta_2=0.999,
        epsilon=1e-07,
        amsgrad=False,
        name="Adam")

    model.compile(
        optimizer=opt,
        loss='binary_crossentropy',
        metrics=['accuracy']
        )
    return model


def plot_history(history, size):
    """
    Plot history of the model training,
    accuracy and loss of the training and validation set
    """
    
    acc = history.history['accuracy']
    val_acc = history.history['val_accuracy']
    loss = history.history['loss']
    val_loss = history.history['val_loss']

    epochs = range(1, len(acc) + 1)

    plt.figure(figsize=(8, 6), dpi=80)

    plt.plot(epochs, acc, 'bo', label='Training acc')
    plt.plot(epochs, val_acc, 'b', label='Validation acc')
    plt.title('Accuracy')
    plt.legend()
    plt.savefig(f"training_acc_{size}.jpg")
    
    plt.figure()

    plt.plot(epochs, loss, 'bo', label='Training loss')
    plt.plot(epochs, val_loss, 'b', label='Validation loss')
    plt.title('Loss')
    plt.legend()
    plt.savefig(f"training_loss_{size}.jpg")


class DataGenerator(Sequence):
    def __init__(self, data_path, labels_path, dataset_size, batch_size, validation_split=0.1, is_validation=False, shuffle=True):
        # preload the encoded numpy data
        # the dataset size is needed to properly load the numpy files
        self.size = dataset_size
            
        self.data = np.memmap(data_path, dtype='float32', mode='r', shape=(self.size, 50, 20, 1))
        self.labels = np.memmap(labels_path, dtype='float32', mode='r', shape=(self.size,))
        self.batch_size = batch_size
        self.shuffle = shuffle
        
        # Determine number of train and validation samples
        self.validation_split = validation_split
        self.num_samples = len(self.data)
        self.num_validation_samples = int(self.num_samples * validation_split)
        self.num_train_samples = self.num_samples - self.num_validation_samples
        
        # Determine indices for validation and training
        indices = np.arange(self.num_samples)
        if shuffle:
            np.random.shuffle(indices)
        
        if is_validation:
            self.indices = indices[self.num_train_samples:]
        else:
            self.indices = indices[:self.num_train_samples]
        
        # Shuffle the data initially
        self.on_epoch_end()

    def __len__(self):
        # Denotes the number of batches per epoch
        return int(np.ceil(len(self.indices) / float(self.batch_size)))

    def __getitem__(self, idx):
        # Generate one batch of data
        batch_indices = self.indices[idx * self.batch_size:(idx + 1) * self.batch_size]
        batch_data = self.data[batch_indices]
        batch_labels = self.labels[batch_indices]
        return batch_data, batch_labels

    def on_epoch_end(self):
        # Updates indices after each epoch for shuffling
        if self.shuffle:
            np.random.shuffle(self.indices)


def train_model(data, labels, dataset_size, ratio, model_file, debug=False):

    # set random state for reproducibility
    random.seed(42)
    np.random.seed(42)
    tf.random.set_seed(42)
    os.environ['TF_DETERMINISTIC_OPS'] = '1'

    train_data_gen = DataGenerator(data, labels, dataset_size, batch_size=32, validation_split=0.1, is_validation=False)
    val_data_gen = DataGenerator(data, labels, dataset_size, batch_size=32, validation_split=0.1, is_validation=True)

    model = compile_model()
    model_history = model.fit(
        train_data_gen,
        validation_data=val_data_gen,
        epochs=10,
        class_weight={0: 1, 1: ratio}
    )

    if debug:
        plot_history(model_history, dataset_size)

    model.save(model_file)


def main():
    parser = argparse.ArgumentParser(description="Train CNN model on encoded miRNA x target binding matrix dataset")
    parser.add_argument('--ratio', type=int, required=True, help="Ratio of pos:neg in the training dataset")
    parser.add_argument('--data', type=str, required=True, help="File with the encoded dataset")
    parser.add_argument('--labels', type=str, required=True, help="File with the dataset labels")
    parser.add_argument('--dataset_size', type=int, required=True, help="Number of samples in the dataset. Needed to properly load the numpy files.")
    parser.add_argument('--model', type=str, required=False, help="Filename to save the trained model")
    parser.add_argument('--debug', type=bool, default=False, help="Set to True to output some plots about training")
    args = parser.parse_args()

    if args.model is None:
        args.model = f"model.keras"

    start = time.time()
    train_model(args.data, args.labels, args.dataset_size, args.ratio, args.model, args.debug)
    end = time.time()
    
    print("Elapsed time: ", end - start, " s.")

if __name__ == "__main__":
    main()


