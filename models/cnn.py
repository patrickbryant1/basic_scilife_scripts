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
from tensorflow.keras.backend import epsilon, clip, sum, log, pow, mean, get_value, set_value, transpose
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

def create_features(df, min_val, max_val):
    '''Get features
    '''
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
    enc_feature1 = []
    enc_feature2 = []

    #Get hmms
    for hgroup in groups:
        group_data = df[df['H_group_x']==hgroup]
        uid1 = [*group_data['uid1']]
        uid2 = [*group_data['uid2']]


        hgroup_s = hgroup.split('.')
        group_name = 'h_'+hgroup_s[0]+'_'+hgroup_s[1]+'_'+hgroup_s[2]+'_'+hgroup_s[3]
        for i in range(0,len(uid1)):
            uids = uid1[i]+'_'+uid2[i]
            enc1_i = pad_cut(enc1[i], 300, 22)
            enc2_i = pad_cut(enc2[i], 300, 22)
            dist = np.asarray([evdist[i]]*22) #Dont want to lose this due to conv
            dist = np.expand_dims(dist, axis=0)
            enc1_i = np.append(enc1_i, dist, axis = 0)
            enc2_i = np.append(enc2_i, dist, axis = 0)


            enc_feature1.append(enc1_i) #Append to list.
            enc_feature2.append(enc2_i) #Append to list.

    #Get RMSDs
    #rmsds = df['RMSD_x'] #rmsds/max_rmsd
    #Bin the TMscore RMSDs
    rmsds = df["RMSD_x"]
    bins = np.arange(min_val, max_val, 0.1)
    binned_rmsd = np.digitize(rmsds, bins)

    X = [np.asarray(enc_feature1),np.asarray(enc_feature2)]
    y = np.eye(45)[binned_rmsd-1] #-1 to start at 0 : Keras needs this (uses indexing)
    #y_binned = np.asarray(binned_rmsds)
    #y_binned = y_binned-1 #Needs to start at 0 for keras
    #y_hot = np.eye(len(bins))[y_binned]

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
max_val = max(complete_df["RMSD_x"])
min_val = min(complete_df["RMSD_x"])
X_train,y_train = create_features(train_df, min_val, max_val)
#Take 500 first points in train
X_train_500 = [[],[]]
y_train_500 = []
count_500 = np.zeros(45)
for i in range(0, len(y_train)):
    pos = np.argmax(y_train[i])
    if count_500[pos] <= 500:
        X_train_500[0].append(X_train[0][i])
        X_train_500[1].append(X_train[1][i])
        y_train_500.append(y_train[i])
        count_500[pos]+=1

X_train_500 = [np.asarray(X_train_500[0]), np.asarray(X_train_500[1])] #convert to arrays
y_train_500 = np.asarray(y_train_500)
X_valid,y_valid = create_features(valid_df, min_val, max_val)
#X_valid = X_valid.reshape(len(X_valid),301,40,1)

#MODEL PARAMETERS
base_epochs = 20
finish_epochs = 2
batch_size = 10
input_dim = X_train_500[0][0].shape
num_classes = max(y_train[0].shape)
seq_length = 301
kernel_size = 6 #they usd 6 and 10 in this paper: https://arxiv.org/pdf/1706.01010.pdf - should do different and cat
filters = 100
drop_rate = 0.5
num_nodes = 300
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
in_1 = keras.Input(shape = input_dim)
in_2 = keras.Input(shape = input_dim)
def resnet(x, num_res_blocks):
	"""Builds a resnet with 1D convolutions of the defined depth.
	"""


    	# Instantiate the stack of residual units
    	#Similar to ProtCNN, but they used batch_size = 64, 2000 filters and kernel size of 21
	for res_block in range(num_res_blocks):
		batch_out1 = BatchNormalization()(x) #Bacth normalize, focus on segment
		activation1 = Activation('relu')(batch_out1)
		conv_out1 = Conv1D(filters = filters, kernel_size = kernel_size, dilation_rate = dilation_rate, input_shape=input_dim, padding ="same")(activation1)
		batch_out2 = BatchNormalization()(conv_out1) #Bacth normalize, focus on segment
		activation2 = Activation('relu')(batch_out2)
		conv_out2 = Conv1D(filters = filters, kernel_size = kernel_size, dilation_rate = dilation_rate, input_shape=input_dim, padding ="same")(activation2)
		x = add([x, conv_out2]) #Skip connection



	return x

#Output (batch, steps(len), filters), filters = channels in next
x1 = resnet(in_1, num_res_blocks)
x2 = resnet(in_2, num_res_blocks)

#Average pool along sequence axis
#x = AveragePooling1D(data_format='channels_first')(x) #data_format='channels_first'
maxpool1 = Lambda(lambda x: keras.backend.max(x1, axis=2))(x1)
#Dense final layer for classification
probabilities = Dense(num_classes, activation='softmax')(maxpool1)

#Model: define inputs and outputs
model = Model(inputs = [in_1, in_2], outputs = probabilities)
sgd = optimizers.SGD(clipnorm=1.)
model.compile(loss='categorical_crossentropy',
              optimizer=sgd,
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

#Fit model
model.fit(X_train_500, y_train_500, batch_size = batch_size,
             epochs=num_epochs,
             validation_data = [X_valid, y_valid],
             shuffle=True, #Dont feed continuously
             callbacks=callbacks)

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
