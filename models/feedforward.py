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
from tensorflow.keras.layers import Dense, Dropout, Activation, Conv1D, Reshape, MaxPooling1D, GlobalAveragePooling1D
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

def create_features(df, h5_path, max_rmsd):
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
    #bins = np.arange(0,4.5,0.1)
    #Bin the TMscore RMSDs
    #binned_rmsds = np.digitize(rmsds, bins)

    #Data
    X = np.asarray(enc_feature)
    y = np.asarray(rmsds/max_rmsd)
    #y_binned = np.asarray(binned_rmsds)
    #y_binned = y_binned-1 #Needs to start at 0 for keras
    #y_hot = np.eye(len(bins))[y_binned]
    #Close h5 file
    h5.close()
    return(X, y)


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

#Max rmsd for normalization
max_rmsd = max(complete_df['RMSD_x'])
X_train,y_train = create_features(train_df, h5_path, max_rmsd)
#X_train = X_train.reshape(len(X_train),301,40,1) #Need 3 dim for 2d conv
X_valid,y_valid = create_features(valid_df, h5_path, max_rmsd)
#X_valid = X_valid.reshape(len(X_valid),301,40,1)

#MODEL PARAMETERS
num_nodes = 300
num_features = 40
input_dim = X_train[0].shape
num_epochs = 2
batch_size = 5

seq_length = 301
kernel_size = 6 #they usd 6 and 10 in paper: https://arxiv.org/pdf/1706.01010.pdf
#MODEL
model = Sequential()
#Should do a resnet with convolutions of 40xsome_number - capture
model.add(Conv1D(seq_length, kernel_size, activation='relu', input_shape=(seq_length, num_features))) #1 filter per residue pair: keras.layers.Conv1D(filters, kernel_size,
model.add(Conv1D(seq_length, kernel_size, activation='relu'))
model.add(MaxPooling1D(3)) #Will be seq_length/pool size outp (they selected 30 numbers in paper --> pool_size = 10)
model.add(Conv1D(128, 3, activation='relu'))
model.add(Conv1D(128, 3, activation='relu'))
model.add(GlobalAveragePooling1D())
model.add(Dropout(0.5))
model.add(Dense(1, kernel_initializer='normal')) #Dense implements the operation: output = activation(dot(input, kernel) + bias)


print(model.summary())
model.compile(loss='mean_absolute_error', #[categorical_focal_loss(alpha=.25, gamma=2)],
              optimizer='adam',
              metrics=['accuracy'])

 #Fit model
model.fit(X_train, y_train, batch_size = batch_size,
             epochs=num_epochs,
             validation_data = [X_valid, y_valid],
             shuffle=True) #Dont feed continuously

pdb.set_trace()
pred = model.predict(X_valid)
average_error = np.average(np.absolute(pred-y_valid))
print(average_error)
#Close h5
h5.close()
pdb.set_trace()
