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
import tables

from tensorflow.keras import regularizers
from tensorflow.keras.constraints import max_norm
from tensorflow.keras.models import Sequential, Model
from tensorflow.keras.layers import Dense, Dropout, Activation, Conv2D, Reshape
from tensorflow.keras.callbacks import TensorBoard

from model_inputs import split_on_h_group, pad_cut
import pdb
#Arguments for argparse module:
parser = argparse.ArgumentParser(description = '''A Neural Network for predicting
                                                RMSD between structural alignments .''')

parser.add_argument('dataframe', nargs=1, type= str,
                  default=sys.stdin, help = 'Path to dataframe in .csv.')

parser.add_argument('h5_path', nargs=1, type= str,
                  default=sys.stdin, help = 'Path to .h5 file with profiles.')

parser.add_argument('out_dir', nargs=1, type= str,
                  default=sys.stdin, help = 'Path to output directory. Include /in end')

#FUNCTIONS
def pad_cut(ar, x):
    '''Pads or cuts a 1D array to len x
    '''
    shape = ar.shape

    if max(shape) < x:
        empty = np.zeros((x,20))
        empty[0:len(ar)]=ar
        ar = empty
    else:
        ar = ar[0:x]

    return ar

def create_features(df, h5_path):
    '''Get features
    '''
    #Open h5
    h5 = tables.open_file(h5_path)
    #Get H_groups
    groups = [*Counter(df['H_group_x']).keys()]
    #Get MLAAdist
    evdist = np.asarray(df['MLAAdist_x'])


    #Save features
    enc_feature = []

    #Get hmms
    for hgroup in groups:
        group_data = df[df['H_group_x']==hgroup]
        uid1 = [*group_data['uid1']]
        uid2 = [*group_data['uid2']]


        hgroup_s = hgroup.split('.')
        group_name = 'h_'+hgroup_s[0]+'_'+hgroup_s[1]+'_'+hgroup_s[2]+'_'+hgroup_s[3]
        for i in range(0,len(uid1)):
            uids = uid1[i]+'_'+uid2[i]
            hmm1 = pad_cut(h5.root[group_name]['hmm1_'+uids][:], 300)
            hmm2 = pad_cut(h5.root[group_name]['hmm2_'+uids][:], 300)
            #tf1 = pad_cut(np.concatenate(h5.root[group_name]['tf1_'+uids][:]), 300*7)
            #tf2 = pad_cut(np.concatenate(h5.root[group_name]['tf2_'+uids][:]), 300*7)
            #ld1 = pad_cut(np.concatenate(h5.root[group_name]['ld1_'+uids][:]), 300*3)
            #ld2 = pad_cut(np.concatenate(h5.root[group_name]['ld2_'+uids][:]), 300*3)

            cat = np.concatenate((hmm1,hmm2), axis = 1)
            dist = np.asarray([evdist[i]]*40)
            dist = np.expand_dims(dist, axis=0)
            cat = np.append(cat, dist, axis = 0)
            enc_feature.append(cat) #Append to list

    #Get RMSDs - should probably normalize with value 4,5 ?
    rmsds = df['RMSD_x']
    bins = np.arange(0,4.5,0.1)
    #Bin the TMscore RMSDs
    binned_rmsds = np.digitize(rmsds, bins)

    #Data
    X = np.asarray(enc_feature)

    y_binned = np.asarray(binned_rmsds)
    y_binned = y_binned-1 #Needs to start at 0 for keras
    y_hot = np.eye(len(bins))[y_binned]
    #Close h5 file
    h5.close()
    return(X, y_hot)


#MAIN
args = parser.parse_args()
dataframe = args.dataframe[0]
h5_path = args.h5_path[0]
out_dir = args.out_dir[0]

#Open h5
h5 = tables.open_file(h5_path)
#Assign data and labels
#Read df
complete_df = pd.read_csv(dataframe)
#Split
train_groups, valid_groups, test_groups = split_on_h_group(complete_df, 0.8)
train_df = complete_df[complete_df['H_group_x'].isin(train_groups)]
valid_df = complete_df[complete_df['H_group_x'].isin(valid_groups)]
test_df = complete_df[complete_df['H_group_x'].isin(test_groups)]

X_train,y_train = create_features(train_df, h5_path)
X_train = X_train.reshape(len(X_train),9000,2,1)
X_valid,y_valid = create_features(valid_df, h5_path)
X_valid = X_valid.reshape(len(X_valid),9000,2,1)
#MODEL PARAMETERS
num_nodes = 300
input_dim = X_train[0].shape
num_labels = max(y_train[0].shape)
num_epochs = 2
batch_size = 20
pdb.set_trace()
#MODEL
model = Sequential()
model.add(Conv2D(64, kernel_size=3, activation='relu', input_shape = input_dim)) #3x3 filter matrix
model.add(Dense(num_nodes, activation="relu"))
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
#Close h5
h5.close()
pdb.set_trace()
