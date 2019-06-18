#! /usr/bin/env python3
# -*- coding: utf-8 -*-


import argparse
import sys
import numpy as np
from ast import literal_eval
import pandas as pd
import glob
from sklearn.model_selection import train_test_split
from sklearn.metrics import classification_report, confusion_matrix
from collections import Counter
import math
import time


from tensorflow.keras import regularizers
from tensorflow.keras.constraints import max_norm
from tensorflow.keras.models import Sequential, Model
from tensorflow.keras.layers import Dense, Dropout, Activation
from tensorflow.keras.callbacks import TensorBoard

import pdb
#Arguments for argparse module:
parser = argparse.ArgumentParser(description = '''A Recurrent Neural Network for predicting
                                                MSD between structural alignments .''')

parser.add_argument('dataframe', nargs=1, type= str,
                  default=sys.stdin, help = 'Path to dataframe in .csv.')

parser.add_argument('out_dir', nargs=1, type= str,
                  default=sys.stdin, help = 'Path to output directory. Include /in end')

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

#MAIN
args = parser.parse_args()
dataframe = args.dataframe[0]
out_dir = args.out_dir[0]

#Assign data and labels
#Read df
complete_df = pd.read_csv(dataframe)

#Get MLAAdist
evdist = complete_df['MLAAdist_x']
evdist = np.asarray(evdist).reshape(len(evdist),1)
#Get encodings
enc1 = []
enc2 = []
[enc1.append(literal_eval(x)) for x in complete_df['enc1']]
[enc2.append(literal_eval(x)) for x in complete_df['enc2']]
#Get longest alignment
enc_lens = []
[enc_lens.append(len(x)) for x in enc1]
#PAD encodings
padlen = max(enc_lens)
enc1 = pad_data(enc1, padlen)
enc2 = pad_data(enc2, padlen)

#Get lengths and reshape to 2D
l1 = complete_df['l1']
l2 = complete_df['l2']
aln_len = complete_df['aln_len']
l1 = np.asarray(l1).reshape(len(l1),1)
l2 = np.asarray(l2).reshape(len(l2),1)
aln_len = np.asarray(aln_len).reshape(len(aln_len),1)

#Create input features
#Concat all features
enc_feature = np.concatenate((np.asarray(enc1), np.asarray(enc2), l1, l2, aln_len, evdist), axis = 1)


#Get RMSDs
rmsds = complete_df['RMSD_x']
bins = np.arange(0,4.5,0.1)
#bins = np.arange(0.5,2.5,0.05)
#bins = np.insert(bins,0, 0)
#bins = np.append(bins, 4.5)
#Bin the TMscore RMSDs
binned_rmsds = np.digitize(rmsds, bins)
#Data
X = np.asarray(enc_feature)
y = np.eye(len(bins))[binned_rmsds-1] #deviations_hot (-1 since start should be 0)

#Split
X_train, X_valid, y_train, y_valid = train_test_split(X, y, test_size=0.2, random_state=42)


#MODEL PARAMETERS
num_nodes = 300
input_dim = len(X_train[0])
num_labels = len(y_train[0])
num_epochs = 20
batch_size = 20
#MODEL
model = Sequential()
model.add(Dense(num_nodes, input_dim=input_dim, activation="relu"))
model.add(Dense(num_nodes/2, activation="relu", kernel_initializer="uniform"))
model.add(Dense(num_labels))
model.add(Activation("softmax"))
model.compile(loss='categorical_crossentropy', #[categorical_focal_loss(alpha=.25, gamma=2)],
              optimizer='adam',
              metrics=['accuracy'])

 #Fit model
model.fit(X_train, y_train, batch_size = batch_size,
             epochs=num_epochs,
             validation_data = [X_valid, y_valid],
             shuffle=True) #Dont feed continuously


pred = np.argmax(model.predict(X_valid), axis = 1)
labels = np.argmax(y_valid, axis = 1)
average_error = np.average(np.absolute(pred-labels))
print(average_error)
pdb.set_trace()
