#!/usr/share/python3
# -*- coding: utf-8 -*-


import argparse
import sys
import numpy as np
from ast import literal_eval
import pandas as pd
import glob
from os import makedirs
from os.path import exists, join

#Preprocessing
from collections import Counter
import math
import time
from ast import literal_eval
from sklearn.model_selection import train_test_split, StratifiedShuffleSplit

#Keras
from tensorflow.keras.models import model_from_json
import tensorflow as tf
from tensorflow.keras import regularizers,optimizers
import tensorflow.keras as keras
from tensorflow.keras.constraints import max_norm
from tensorflow.keras.callbacks import ModelCheckpoint, LearningRateScheduler, Callback
from tensorflow.keras.models import Sequential, Model
from tensorflow.keras.layers import Dense, Dropout, Activation, Conv1D, Reshape, MaxPooling1D, Dot, Masking
from tensorflow.keras.layers import Activation, RepeatVector, Permute, Multiply, Lambda, GlobalAveragePooling1D
from tensorflow.keras.layers import concatenate, add, Conv1D, BatchNormalization, Flatten, Subtract
from tensorflow.keras.backend import epsilon, clip, get_value, set_value, transpose, variable, square
from tensorflow.layers import AveragePooling1D
from tensorflow.keras.losses import mean_absolute_error, mean_squared_error

import pdb


#Arguments for argparse module:
parser = argparse.ArgumentParser(description = '''A softmax classifier that uses embeddings generated from a pretrained model.''')

parser.add_argument('json_file', nargs=1, type= str,
                  default=sys.stdin, help = 'path to .json file with keras model to be opened')

parser.add_argument('weights', nargs=1, type= str,
                  default=sys.stdin, help = '''path to .h5 file containing weights for net.''')

parser.add_argument('encodings', nargs=1, type= str,
                  default=sys.stdin, help = 'Path to np array with encoded aa sequences.')

parser.add_argument('dataframe', nargs=1, type= str,
                  default=sys.stdin, help = '''path to data to be used for prediction.''')

parser.add_argument('class_embeddings', nargs=1, type= str,
                  default=sys.stdin, help = '''path to class embeddings.''')

parser.add_argument('out_dir', nargs=1, type= str,
                  default=sys.stdin, help = '''path to output directory.''')




#FUNCTIONS

def pad_cut(ar, x, y):
    '''Pads or cuts a 1D array to len x
    '''
    shape = ar.shape

    if max(shape) < x:
        empty = np.zeros((x,y))
        empty[0:len(ar)]=ar
        ar = empty
    else:
        ar = ar[0:x]
    return ar

def load_model(json_file, weights):

	global model

	json_file = open(json_file, 'r')
	model_json = json_file.read()
	model = model_from_json(model_json)
	model.load_weights(weights)
	model._make_predict_function()
	return model

#MAIN
args = parser.parse_args()
json_file = (args.json_file[0])
weights = (args.weights[0])
encodings = args.encodings[0]
dataframe = args.dataframe[0]
class_embeddings = np.load(args.class_embeddings[0], allow_pickle=True)
out_dir = args.out_dir[0]

#Assign data and labels
#Read df
df = pd.read_csv(dataframe)

#Assign data and labels
#Read data
X = np.load(encodings, allow_pickle=True)
y = np.asarray([*df['group_enc']])
y_hot = np.eye(max(y)+1)[y]
#Pad X
padded_X = []
for i in range(0,len(X)):
    padded_X.append(pad_cut(X[i], 300, 21))
X = np.asarray(padded_X)


#Load and run model
model = load_model(json_file, weights)

#Get embedding layers
features1 = Model(inputs=model.input, outputs=model.get_layer('features1').output)
features2 = Model(inputs=model.input, outputs=model.get_layer('features2').output)

#Get average embeddings for all entries
emb1 = np.asarray(features1.predict([X,X]))
emb2 = np.asarray(features2.predict([X,X]))
average_emb = np.average([emb1, emb2], axis = 0)
pdb.set_trace()
#Compute L1 distance to all class_embeddings
L_dists = []
for i in range(0, len(average_emb)):
    true = y[i]
    diff = np.absolute(class_embeddings-average_emb[i])
    L_dists.append(diff)
pdb.set_trace()

#MODEL
x = keras.Input(shape = input_dim)
#Create softmax classifier
probability = Dense(1, activation='softmax')(x)



softmax_clf = Model(inputs = [x], outputs = probability)
opt = optimizers.Adam(clipnorm=1.)
softmax_clf.compile(loss='categorical_crossentropy',
              metrics = ['accuracy'],
              optimizer=opt)
softmax_clf.fit(L_dists, y_hot, batch_size = 1)
pdb.set_trace()
