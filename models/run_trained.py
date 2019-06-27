#! /usr/bin/env python3
# -*- coding: utf-8 -*-


import argparse
import sys
import numpy as np
import pandas as pd
from ast import literal_eval
from sklearn.model_selection import train_test_split
from tensorflow.keras.models import model_from_json
import math
import time
import tensorflow as tf
import tensorflow.keras as keras
from tensorflow.keras import regularizers, backend
from tensorflow.keras.constraints import max_norm
from tensorflow.keras.models import Sequential, Model
from tensorflow.keras.layers import Dense, Embedding, Flatten
from tensorflow.keras.layers import Bidirectional,CuDNNLSTM, Dropout, BatchNormalization
from tensorflow.keras.layers import Reshape, Activation, RepeatVector, Permute, multiply, Lambda
from tensorflow.keras.layers import concatenate, add, Conv1D
from tensorflow.keras.callbacks import ModelCheckpoint, LearningRateScheduler, Callback
from tensorflow.keras.preprocessing.sequence import TimeseriesGenerator
from tensorflow.keras.callbacks import TensorBoard

import pdb

from model_inputs import split_on_h_group
#Arguments for argparse module:
parser = argparse.ArgumentParser(description = '''A program that reads a keras model from a .json and a .h5 file''')

parser.add_argument('json_file', nargs=1, type= str,
                  default=sys.stdin, help = 'path to .json file with keras model to be opened')

parser.add_argument('weights', nargs=1, type= str,
                  default=sys.stdin, help = '''path to .h5 file containing weights for net.''')

parser.add_argument('dataframe', nargs=1, type= str,
                  default=sys.stdin, help = '''path to data to be used for prediction.''')
#FUNCTIONS

def pad_data(X, padlen):
	'''Pads entries in each batch with zeros to have equal lengths
	'''

	#Loop through X
	X_pad = [] #save padded data
	for i in range(0,len(X)):
		if len(X[i])>padlen:
			pdb.set_trace()
		X_pad.append(np.pad(X[i], (0,padlen-len(X[i])), 'constant'))

	return X_pad

def load_model(json_file, weights):

	global model

	json_file = open(json_file, 'r')
	model_json = json_file.read()
	model = model_from_json(model_json)
	model.load_weights(weights)
	model._make_predict_function()
	model.compile(optimizer='adam', loss='categorical_crossentropy', metrics=['accuracy'])
	return model

#MAIN
args = parser.parse_args()
json_file = (args.json_file[0])
weights = (args.weights[0])
dataframe = args.dataframe[0]

#Assign data and labels
#Read df
complete_df = pd.read_csv(dataframe)
#Get encodings
enc1 = []
enc2 = []
[enc1.append(literal_eval(x)) for x in complete_df['enc1']]
[enc2.append(literal_eval(x)) for x in complete_df['enc2']]
#Get longest alignment
enc_lens = []
[enc_lens.append(len(x)) for x in enc1]
#Get MLAAdist
evdist = complete_df['MLAAdist_x']
#Convert to array
X = [np.asarray(enc1), np.asarray(enc2), np.asarray(evdist)]
X = np.asarray(X)#Convert to np array
X = X.T #transpose
#One-hot encode binned data
#-4.2478745796221045 2.871030140344752
bins = np.arange(0.5,2.5,0.25)
bins = np.insert(bins,0, 0)
bins = np.append(bins, 4.5)
#Bin the TMscore RMSDs
rmsds = complete_df['RMSD_x']


binned_rmsds = np.digitize(rmsds, bins)
np.asarray(binned_rmsds)
y = np.eye(len(bins))[binned_rmsds] #deviations_hot

#Split train data to use 80% for training and 10% for validation and 10 % for testing.
X_train, X_valid, y_train, y_valid = train_test_split(X, y, test_size=0.2, random_state=42)
#Random state = 42 guarantees the split is the same every time. This can be both bad and good, depending on
#the selction. It makes the split consistent across changes to the network though.
#Get test data
X_valid, X_test, y_valid, y_test = train_test_split(X_valid, y_valid, test_size=0.5, random_state=42)


#Transpose back
X_train = X_train.T
X_valid = X_valid.T
X_test = X_test.T

#Pad data
#Get longest alignment
padlen = max(enc_lens)
#Set as two separate arrays for encoding representation
X_train_1 = pad_data(X_train[0], padlen)
X_train_2 = pad_data(X_train[1], padlen)
X_train = [np.asarray(X_train_1),np.asarray(X_train_2)] #, np.asarray(X_train_3)]
X_valid_1 = pad_data(X_valid[0], padlen)
X_valid_2 = pad_data(X_valid[1], padlen)
X_valid = [np.asarray(X_valid_1),np.asarray(X_valid_2)]

#Load and run model
model = load_model(json_file, weights)
pred = model.predict(X_valid)

argmax_pred = tf.argmax(pred, 1)

sess = tf.Session()
called = sess.run(argmax_pred)

argmax_valid = np.argmax(y_valid, axis = 1)
average_error = np.average(np.absolute(called-argmax_valid))
print(average_error)
pdb.set_trace()
