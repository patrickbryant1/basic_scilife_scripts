#! /usr/bin/env python3
# -*- coding: utf-8 -*-


import argparse
import sys
import numpy as np
from ast import literal_eval
import pandas as pd
import glob

#Preprocessing
from sklearn.metrics import classification_report, confusion_matrix
from collections import Counter
import math
import time
import tables
from ast import literal_eval

#Keras
from tensorflow.keras import regularizers
import tensorflow.keras as keras
from tensorflow.keras.constraints import max_norm
from tensorflow.keras.models import Sequential, Model
from tensorflow.keras.layers import Dense, Dropout, Activation, Conv1D, Reshape, MaxPooling1D
from tensorflow.keras.layers import Activation, RepeatVector, Permute, multiply, Lambda, GlobalAveragePooling1D
from tensorflow.keras.layers import concatenate, add, Conv1D, BatchNormalization, Flatten
from tensorflow.keras.callbacks import TensorBoard

#visualization
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

def create_features(df, h5_path, max_rmsd):
    '''Get features
    '''
    #Open h5
    h5 = tables.open_file(h5_path)
    #Get H_groups
    groups = [*Counter(df['H_group_x']).keys()]
    #Get MLAAdist
    evdist = np.asarray(df['MLAAdist_x'])

    #Get encodings
    enc1 = []
    enc2 = []
    [enc1.append(np.eye(22)[literal_eval(x)]) for x in [*df['enc1']]]
    [enc2.append(np.eye(22)[literal_eval(x)]) for x in [*df['enc2']]]
    #onehot
    #enc1_hot = np.eye(len(enc1[0]))[enc1]
    #enc2_hot = np.eye(len(enc1[0]))[enc2]

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
            hmm1 = pad_cut(h5.root[group_name]['hmm1_'+uids][:], 300, 20)
            hmm2 = pad_cut(h5.root[group_name]['hmm2_'+uids][:], 300, 20)
            #tf1 = pad_cut(np.concatenate(h5.root[group_name]['tf1_'+uids][:]), 300*7)
            #tf2 = pad_cut(np.concatenate(h5.root[group_name]['tf2_'+uids][:]), 300*7)
            #ld1 = pad_cut(np.concatenate(h5.root[group_name]['ld1_'+uids][:]), 300*3)
            #ld2 = pad_cut(np.concatenate(h5.root[group_name]['ld2_'+uids][:]), 300*3)
            enc1_i = pad_cut(enc1[i], 300, 22)
            enc2_i = pad_cut(enc2[i], 300, 22)
            cat = np.concatenate((hmm1,hmm2, enc1_i, enc2_i), axis = 1)
            dist = np.asarray([evdist[i]]*84)
            dist = np.expand_dims(dist, axis=0)
            cat = np.append(cat, dist, axis = 0)
            enc_feature.append(cat) #Append to list


    #Get RMSDs
    #rmsds = df['RMSD_x'] #rmsds/max_rmsd
    #bins = np.arange(0,4.5,0.1)
    #Bin the TMscore RMSDs
    #binned_rmsds = np.digitize(rmsds, bins)
    deviations = [*df['deviation']]
    bins = np.arange(min(deviations),max(deviations),0.2)
    binned_deviations = np.digitize(deviations, bins)
    #t= 0.15 #Maybe I should bin the deviations more than just +/- t: likely those deviating more are more different
    #for i in range(0, len(deviations)):
    #    if deviations[i] > t:
    #        binned_deviations.append(2)
    #    if deviations[i] < -t:
    #        binned_deviations.append(0)
    #    if np.absolute(deviations[i]) <= t:
    #        binned_deviations.append(1)
    deviations_hot = np.eye(3)[binned_deviations]


    #Data
    X = np.asarray(enc_feature)
    y = deviations_hot
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
num_features = 84 #Perhaps add a one if not gap for each reisude = 42 features
input_dim = X_train[0].shape
num_epochs = 10
batch_size = 5
num_classes = 3
seq_length = 301
kernel_size = 1 #they usd 6 and 10 in paper: https://arxiv.org/pdf/1706.01010.pdf
drop_rate = 0.5
num_nodes = 301

#MODEL
in_params = keras.Input(shape = input_dim)

def resnet(x, num_res_blocks, num_classes=num_classes):
	"""Builds a resnet with 1D convolutions of the defined depth.
	"""

# Instantiate the stack of residual units
	for res_block in range(num_res_blocks):
		conv_out1 = Conv1D(seq_length, kernel_size, activation='relu', input_shape=(seq_length, num_features))(x)
		conv_out1 = BatchNormalization()(conv_out1) #Bacth normalize, focus on segment
		conv_out1 = Dropout(rate = drop_rate)(conv_out1) #Dropout
		conv_out2 = Conv1D(seq_length, kernel_size, activation='relu', input_shape=(seq_length, num_features))(conv_out1)
		conv_out2 = BatchNormalization()(conv_out2) #Bacth normalize, focus on segment
		conv_out2 = Dropout(rate = drop_rate)(conv_out2) #Dropout

		conv_add1 = add([conv_out1, conv_out2]) #Skip connection, add before or after dropout?

		conv_out3 = Conv1D(seq_length, kernel_size, activation='relu', input_shape=(seq_length, num_features))(conv_add1)
		conv_out3 = BatchNormalization()(conv_out3) #Bacth normalize, focus on segment
		conv_out3 = Dropout(rate = drop_rate)(conv_out3) #Dropout
		x = add([conv_out2, conv_out3])

	return x
#Create resnet and get outputs
x = resnet(in_params, 1, num_classes)


#Attention layer - information will be redistributed in the backwards pass
attention = Dense(1, activation='tanh')(x) #Normalize and extract info with tanh activated weight matrix (hidden attention weights)
attention = Flatten()(attention) #Make 1D
attention = Activation('softmax')(attention) #Softmax on all activations (normalize activations)
attention = RepeatVector(num_nodes)(attention) #Repeats the input "num_nodes" times.
attention = Permute([2, 1])(attention) #Permutes the dimensions of the input according to a given pattern. (permutes pos 2 and 1 of attention)

sent_representation = multiply([x, attention]) #Multiply input to attention with normalized activations
sent_representation = Lambda(lambda xin: keras.backend.sum(xin, axis=-2), output_shape=(num_nodes*2,))(sent_representation) #Sum all attentions

#Dense final layer for classification
probabilities = Dense(num_classes, activation='softmax')(sent_representation)

#Model: define inputs and outputs
model = Model(inputs = in_params, outputs = probabilities)

model.compile(loss='categorical_crossentropy', #[categorical_focal_loss(alpha=.25, gamma=2)],
              optimizer='adam',
              metrics=['accuracy'])
#Summary of model
print(model.summary())

 #Fit model
model.fit(X_train, y_train, batch_size = batch_size,
             epochs=num_epochs,
             validation_data = [X_valid, y_valid],
             shuffle=True) #Dont feed continuously


pred = model.predict(X_valid)
average_error = np.average(np.absolute(pred-y_valid))*max_rmsd
print(average_error)
#Close h5
h5.close()
pdb.set_trace()
