#! /usr/bin/env python3
# -*- coding: utf-8 -*-


import argparse
import sys
import numpy as np
from ast import literal_eval
import pandas as pd
import glob

#Preprocessing
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
from tensorflow.keras.layers import Dense, Dropout, Activation, Conv1D, Reshape, MaxPooling1D, dot, Masking
from tensorflow.keras.layers import Activation, RepeatVector, Permute, Multiply, Lambda, GlobalAveragePooling1D
from tensorflow.keras.layers import concatenate, add, Conv1D, BatchNormalization, Flatten, Subtract
from tensorflow.keras.backend import epsilon, clip, get_value, set_value, transpose, variable, square
from tensorflow.layers import AveragePooling1D
from tensorflow.keras.losses import mean_absolute_error, mean_squared_error
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

parser.add_argument('params_file', nargs=1, type= str,
                  default=sys.stdin, help = 'Path to file with net parameters')

parser.add_argument('out_dir', nargs=1, type= str,
                  default=sys.stdin, help = 'Path to output directory. Include /in end')


#FUNCTIONS
def read_net_params(params_file):
    '''Read and return net parameters
    '''
    net_params = {} #Save information for net

    with open(params_file) as file:
        for line in file:
            line = line.rstrip() #Remove newlines
            line = line.split("=") #Split on "="

            net_params[line[0]] = line[1]


    return net_params

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

def prob_bin(scores, bins):
    '''Bin scores into probabilities in different bins
    '''

    pb_scores = []
    bin_size = bins[1]-bins[0]
    for i in scores:
        pb_score = np.zeros(len(bins))
        if i == 1.0:
            pb_score[-1] = 1

        if float(i) != 0.0 and float(i) != 1.0:
            for j in range(0, len(bins)):
                if i < bins[j+1]:
                    stri = str(i)
                    num_decimals = len(stri)-2
                    n2 = stri[len(str(bin_size)):]
                    pb_score[j] = (1-(float('0.'+'0'*(num_decimals-2)+n2)/bin_size))
                    pb_score[j+1] = (float('0.'+'0'*(num_decimals-2)+n2)/bin_size)
                    break #found which bin

        pb_scores.append(pb_score)
    return np.asarray(pb_scores)

def create_features(df,seq_length, bins):
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
            enc1_i = pad_cut(enc1[i], seq_length, 22)
            enc2_i = pad_cut(enc2[i], seq_length, 22)

            enc_feature1.append(enc1_i) #Append to list.
            enc_feature2.append(enc2_i) #Append to list.

    #Get LDDT
    scores = np.asarray(df["global_lddt"])
    pb_scores = prob_bin(scores, bins)


    X = [np.asarray(enc_feature1),np.asarray(enc_feature2)]
    y = np.asarray(pb_scores)
    return(X, y)


#MAIN
args = parser.parse_args()
dataframe = args.dataframe[0]
params_file = args.params_file[0]
out_dir = args.out_dir[0]

#Assign data and labels
#Read df
complete_df = pd.read_csv(dataframe)
#Split
train_groups, valid_groups, test_groups = split_on_h_group(complete_df, 0.8)
train_df = complete_df[complete_df['H_group_x'].isin(train_groups)]
valid_df = complete_df[complete_df['H_group_x'].isin(valid_groups)]
test_df = complete_df[complete_df['H_group_x'].isin(test_groups)]

#Max rmsd for normalization
max_val = max(complete_df["global_lddt"])
min_val = min(complete_df["global_lddt"])
bins = np.arange(min_val, max_val+0.05, 0.05)


seq_length = 300
#Make the model mix what is fed to x1 and x2 - so both convnets learn the same thing!
X_train,y_train = create_features(train_df, seq_length, bins)

#Take p first points in train
X_train_p = [[],[]]
y_train_p = []
count_p = np.zeros(len(bins))
p = 1000
for i in range(0, len(y_train)):
      pos = np.argmax(y_train[i])
      if count_p[pos] <= p:
          X_train_p[0].append(X_train[0][i])
          X_train_p[1].append(X_train[1][i])
          y_train_p.append(y_train[i])
          count_p[pos]+=1

X_train_p = [np.asarray(X_train_p[0]), np.asarray(X_train_p[1])] #convert to arrays
y_train_p = np.asarray(y_train_p)
#y_train_p = np.matmul(y_train_p, bins)
X_valid,y_valid = create_features(valid_df, seq_length, bins)
#y_valid = np.matmul(y_valid, bins)
bins = np.expand_dims(bins, axis=0)
#Tensorboard for logging and visualization
log_name = str(time.time())
tensorboard = TensorBoard(log_dir=out_dir+log_name)

######MODEL######
#Parameters
#net_params = read_net_params(params_file)
base_epochs = 10
finish_epochs = 3
batch_size = 8
input_dim = X_train[0][0].shape
num_classes = max(bins.shape)
kernel_size = 21 #google uses 21
filters = 100
drop_rate = 0.5
num_nodes = 300
num_res_blocks = 3
dilation_rate = 3

#lr opt
find_lr = False
#LR schedule
num_epochs = base_epochs+finish_epochs
max_lr = 0.001
min_lr = max_lr/10
lr_change = (max_lr-min_lr)/(base_epochs/2-1) #Reduce further last epochs

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
        #Downsample - half filters
		conv_out2 = Conv1D(filters = int(filters/2), kernel_size = kernel_size, dilation_rate = 1, input_shape=input_dim, padding ="same")(activation2)
		x = Conv1D(filters = int(filters/2), kernel_size = kernel_size, dilation_rate = 1, input_shape=input_dim, padding ="same")(x)
		x = add([x, conv_out2]) #Skip connection



	return x

#Initial convolution
in_1_conv = Conv1D(filters = filters, kernel_size = kernel_size, dilation_rate = 2, input_shape=input_dim, padding ="same")(in_1)
in_2_conv = Conv1D(filters = filters, kernel_size = kernel_size, dilation_rate = 2, input_shape=input_dim, padding ="same")(in_2)
#Output (batch, steps(len), filters), filters = channels in next
x1 = resnet(in_1_conv, num_res_blocks)
x2 = resnet(in_2_conv, num_res_blocks)

#Average pool along sequence axis
#x = AveragePooling1D(data_format='channels_first')(x) #data_format='channels_first'
maxpool1 = MaxPooling1D(pool_size=seq_length)(x1)
maxpool2 = MaxPooling1D(pool_size=seq_length)(x2)
#cat = concatenate([maxpool1, maxpool2]) #Cat convolutions

flat1 = Flatten()(maxpool1)  #Flatten
flat2 = Flatten()(maxpool2)  #Flatten

 # Add a customized layer to compute the absolute difference between the encodings
L1_layer = Lambda(lambda tensors:abs(tensors[0] - tensors[1]))
L1_distance = L1_layer([flat1, flat2])

#L2-normalize samples along the dot product axis before taking the dot product.
#cos_dist = keras.layers.Dot(axes = 1, normalize=True)([flat1, flat2]) #normalize = True means the cosine similarity is calculated

#Dense final layer for classification
probabilities = Dense(num_classes, activation='softmax')(L1_distance)

#Custom loss
def bin_loss(y_true, y_pred):
     bins_K = variable(value=bins)
     pred_vals = dot([y_pred, bins_K], axes = 1)
     #true_vals = dot([y_true, bins_K], axes = 1)

     loss = mean_absolute_error(y_true, pred_vals)
     return loss


#Model: define inputs and outputs
model = Model(inputs = [in_1, in_2], outputs = probabilities)
sgd = optimizers.SGD(clipnorm=1.)
model.compile(loss='categorical_crossentropy',
              optimizer=sgd)

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
#Should shuffle uid1 and uid2 in X[0] vs X[1]
model.fit(X_train_p, y_train_p, batch_size = batch_size,
             epochs=num_epochs,
             validation_data = [X_valid, y_valid],
             shuffle=True, #Dont feed continuously
             callbacks=callbacks)

pred = model.predict(X_valid)
pred_rmsds = np.matmul(pred, bins.T)
true = np.matmul(y_valid,bins.T)
mean_error = np.average(np.absolute(pred_rmsds-true))
#Prind validation predictions to file for further analysis
with open('validation.tsv', 'w') as file:
    for i in range(0, len(pred)):
        file.write(str(pred_rmsds[i])+'\t'+str(true[i])+'\n')
pdb.set_trace()
