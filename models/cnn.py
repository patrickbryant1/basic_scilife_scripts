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
from tensorflow.keras import regularizers,optimizers
import tensorflow.keras as keras
from tensorflow.keras.constraints import max_norm
from tensorflow.keras.callbacks import ModelCheckpoint, LearningRateScheduler, Callback
from tensorflow.keras.models import Sequential, Model
from tensorflow.keras.layers import Dense, Dropout, Activation, Conv1D, Reshape, MaxPooling1D
from tensorflow.keras.layers import Activation, RepeatVector, Permute, multiply, Lambda, GlobalAveragePooling1D
from tensorflow.keras.layers import concatenate, add, Conv1D, BatchNormalization, Flatten
from tensorflow.keras.backend import epsilon, clip, sum, log, pow, mean, get_value, set_value
from tensorflow.layers import AveragePooling1D
#visualization
from tensorflow.keras.callbacks import TensorBoard
#Custom
from model_inputs import split_on_h_group, pad_cut
from lr_finder import LRFinder
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

def create_features(df, h5_path, min_val, max_val):
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
            #hmm1 = pad_cut(h5.root[group_name]['hmm1_'+uids][:], 300, 20)
            #hmm2 = pad_cut(h5.root[group_name]['hmm2_'+uids][:], 300, 20)
            #tf1 = pad_cut(np.concatenate(h5.root[group_name]['tf1_'+uids][:]), 300*7)
            #tf2 = pad_cut(np.concatenate(h5.root[group_name]['tf2_'+uids][:]), 300*7)
            #ld1 = pad_cut(np.concatenate(h5.root[group_name]['ld1_'+uids][:]), 300*3)
            #ld2 = pad_cut(np.concatenate(h5.root[group_name]['ld2_'+uids][:]), 300*3)
            enc1_i = pad_cut(enc1[i], 300, 22)
            enc2_i = pad_cut(enc2[i], 300, 22)
            cat = np.concatenate((enc1_i, enc2_i), axis = 1)
            dist = np.asarray([evdist[i]]*min(cat[0].shape))
            dist = np.expand_dims(dist, axis=0)
            cat = np.append(cat, dist, axis = 0)
            enc_feature.append(cat.T) #Append to list


    #Get RMSDs
    #rmsds = df['RMSD_x'] #rmsds/max_rmsd
    #bins = np.arange(0,4.5,0.1)
    #Bin the TMscore RMSDs
    #binned_rmsds = np.digitize(rmsds, bins)
    deviations = [*df['deviation']]
    binned_deviations = []
    t= 0.15 #Maybe I should bin the deviations more than just +/- t: likely those deviating more are more different
    for i in range(0, len(deviations)):
        if deviations[i] > t:
            binned_deviations.append(2)
        if deviations[i] < -t:
            binned_deviations.append(0)
        if np.absolute(deviations[i]) <= t:
            binned_deviations.append(1)
    deviations_hot = np.eye(3)[binned_deviations]
    #Counter({0: 13099, 1: 12521, 2: 11575})
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
max_val = max(complete_df['deviation'])
min_val = min(complete_df['deviation'])
X_train,y_train = create_features(train_df, h5_path, min_val, max_val)
#X_train = X_train.reshape(len(X_train),301,40,1) #Need 3 dim for 2d conv
X_valid,y_valid = create_features(valid_df, h5_path, min_val, max_val)
#X_valid = X_valid.reshape(len(X_valid),301,40,1)

#MODEL PARAMETERS
num_features = min(X_train[0].shape) #Perhaps add a one if not gap for each reisude = 42 features
input_dim = X_train[0].shape
base_epochs = 10
finish_epochs = 2
batch_size = 10
num_classes = max(y_train[0].shape)
seq_length = 301
kernel_size = 1 #they usd 6 and 10 in this paper: https://arxiv.org/pdf/1706.01010.pdf - but then no dilated conv
filters = 200
drop_rate = 0.5
num_nodes = 301
num_res_blocks = 2
dilation_rate = 3

#lr opt
find_lr = False
#LR schedule

num_epochs = base_epochs+finish_epochs
max_lr = 0.001
min_lr = max_lr/10
lr_change = (max_lr-min_lr)/(base_epochs/2-1) #Reduce further lst three epochs

#MODEL
in_params = keras.Input(shape = input_dim)

def resnet(x, num_res_blocks):
	"""Builds a resnet with 1D convolutions of the defined depth.
	"""

# Instantiate the stack of residual units
#Similar to ProtCNN, but they used batch_size = 64, 2000 filters and kernel size of 21
	for res_block in range(num_res_blocks):
		batch_out1 = BatchNormalization()(x) #Bacth normalize, focus on segment
		relu_out1 = Dense(num_nodes, activation='relu')(batch_out1)
        #input_shape=(10, 128) for time series sequences of 10 time steps with 128 features per step
		conv_out1 = Conv1D(filters = filters, kernel_size = kernel_size, activation='relu', dilation_rate = dilation_rate, input_shape=input_dim)(relu_out1)
		#conv_out1 = Dropout(rate = drop_rate)(conv_out1) #Dropout
		batch_out2 = BatchNormalization()(conv_out1) #Bacth normalize, focus on segment
		relu_out2 = Dense(filters, activation='relu')(batch_out2)
        #Filters to match input dim
		conv_out2 = Conv1D(filters = 301, kernel_size = kernel_size, activation='relu', input_shape=input_dim)(relu_out2)


		conv_add1 = add([x, conv_out2]) #Skip connection

		x = conv_add1

	return conv_add1
#Create resnet and get outputs
x = resnet(in_params, num_res_blocks)

#Average pool along sequence axis
#x = AveragePooling1D(data_format='channels_first')(x) #data_format='channels_first'
avgpool = Lambda(lambda x: keras.backend.max(x, axis=2))(x)
#Dense final layer for classification
#flat = Flatten()(x)
probabilities = Dense(num_classes, activation='softmax')(avgpool)

#Model: define inputs and outputs
model = Model(inputs = in_params, outputs = probabilities)
sgd = optimizers.SGD(clipnorm=1.)
model.compile(loss='categorical_crossentropy',
              optimizer='SGD',
              metrics=['accuracy'])

#LR schedule
def lr_schedule(epochs):
  '''lr scheduel according to one-cycle policy.
  '''

  #Increase lrate in beginning
  if epochs == 0:
    lrate = min_lr
  elif (epochs <(base_epochs/2) and epochs > 0):
    lrate = min_lr+(epochs*lr_change)
  #Decrease further below min_lr last three epochs
  elif epochs >= base_epochs:
    lrate = min_lr/(2*(epochs+1-base_epochs))
  #After the max lrate is reached, decrease it back to the min
  else:
    lrate = max_lr-((epochs-(base_epochs/2))*lr_change)

  print(epochs,lrate)
  return lrate


if find_lr == True:
    num_epochs = 1
    batch = 0
    class LossHistory(Callback):
        def on_train_begin(self, logs={}):
            self.losses = []
            self.lrs = []

        def on_batch_end(self, batch, logs={}):
            batch +=1
            self.losses.append(logs.get('loss'))
            set_value(model.optimizer.lr, 0.000001+(batch*((1-0.000001)/(len(X_train)/batch_size))))
            lr = get_value(self.model.optimizer.lr)
            self.lrs.append(lr)

    history = LossHistory()
    callbacks=[history]
else:
  lrate = LearningRateScheduler(lr_schedule)
  callbacks=[lrate]
  validation_data=(X_valid, y_valid)


#Summary of model
print(model.summary())

#Class weights - weight harder classes more
class_weight = {0: 1.,
                1: 1.3,
                2: 1.5}
#Fit model
model.fit(X_train, y_train, batch_size = batch_size,
             epochs=num_epochs,
             validation_data = [X_valid, y_valid],
             shuffle=True, #Dont feed continuously
             callbacks=callbacks,
             class_weight = class_weight)

#Get history: print(history.losses)
if find_lr == True:
    lrs = history.lrs
    losses = history.losses
    with open('lr_plot.tsv', 'w') as file:
        for i in range(0, len(lrs)):
            file.write(str(lrs[i])+'\t'+str(losses[i])+'\n')

pred = np.argmax(model.predict(X_valid), axis = 1)
y_valid = np.argmax(y_valid, axis = 1)
average_error = np.average(np.absolute(pred-y_valid))
print(average_error)
#Prind validation predictions to file
with open('validation.tsv', 'w') as file:
    for i in range(0, len(pred)):
        file.write(str(pred[i])+'\t'+str(y_valid[i])+'\n')
pdb.set_trace()
