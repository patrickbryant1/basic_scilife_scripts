#!/usr/share/python
# -*- coding: utf-8 -*-


import argparse
import sys
import numpy as np
from ast import literal_eval
import pandas as pd
import glob

#Preprocessing
#from sklearn.metrics import classification_report, confusion_matrix
from collections import Counter
import math
import time
#import tables
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
            #dist = np.asarray([evdist[i]]*22) #Dont want to lose this due to conv
            #dist = np.expand_dims(dist, axis=0)
            #enc1_i = np.append(enc1_i, dist, axis = 0)
            #enc2_i = np.append(enc2_i, dist, axis = 0)


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
out_dir = args.out_dir[0]

#Assign data and labels
#Read df
complete_df = pd.read_csv(dataframe)
#Split
train_groups, valid_groups, test_groups = split_on_h_group(complete_df, 0.8)
train_df = complete_df[complete_df['H_group_x'].isin(train_groups)]
valid_df = complete_df[complete_df['H_group_x'].isin(valid_groups)]
test_df = complete_df[complete_df['H_group_x'].isin(test_groups)]


#Tensorboard for logging and visualization
log_name = str(time.time())
tensorboard = TensorBoard(log_dir=out_dir+log_name)

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
base_epochs = 40
finish_epochs = 2
batch_size = 32
input_dim = X_train_500[0][0].shape
num_classes = max(y_train[0].shape)
seq_length = 300
kernel_size = 21 #google uses 21
filters = 1100
drop_rate = 0.5
num_nodes = 300
num_res_blocks = 1
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
		#Bottleneck convolution: downsample to reduce channels
		conv_out2 = Conv1D(filters = int(filters/2), kernel_size = kernel_size, dilation_rate = dilation_rate, input_shape=input_dim, padding ="same")(activation2)
		x = Conv1D(filters = int(filters/2), kernel_size = kernel_size, dilation_rate = 1, input_shape=input_dim, padding ="same")(x)
		x = add([x, conv_out2]) #Skip connection



	return x

#Initial convolution - to make shapes match
in_1_conv = Conv1D(filters = filters, kernel_size = kernel_size, dilation_rate = 2, input_shape=input_dim, padding ="same")(in_1)
in_2_conv = Conv1D(filters = filters, kernel_size = kernel_size, dilation_rate = 2, input_shape=input_dim, padding ="same")(in_2)
#Output (batch, steps(len), filters), filters = channels in next
x1 = resnet(in_1_conv, num_res_blocks)
x2 = resnet(in_2_conv, num_res_blocks)

#Average pool along sequence axis
#x = AveragePooling1D(data_format='channels_first')(x) #data_format='channels_first'
maxpool1 = MaxPooling1D(pool_size=seq_length)(x1)
maxpool2 = MaxPooling1D(pool_size=seq_length)(x2)
cat = concatenate([maxpool1, maxpool2]) #Cat convolutions
flat = Flatten()(cat) #Flatten for dense layer
#Dense final layer for classification
probabilities = Dense(num_classes, activation='softmax')(flat)

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

#Lrate
lrate = LearningRateScheduler(lr_schedule)
callbacks=[lrate, tensorboard]



#Summary of model
print(model.summary())

#Fit model
model.fit(X_train_500, y_train_500, batch_size = batch_size,
             epochs=num_epochs,
             validation_data = [X_valid, y_valid],
             shuffle=True, #Dont feed continuously
             callbacks=callbacks)


pred = np.argmax(model.predict(X_valid), axis = 1)
y_valid = np.argmax(y_valid, axis = 1)
average_error = np.average(np.absolute(pred-y_valid))
print(average_error)
#Prind validation predictions to file
with open(out_dir+'validation.tsv', 'w') as file:
    for i in range(0, len(pred)):
        file.write(str(pred[i])+'\t'+str(y_valid[i])+'\n')
