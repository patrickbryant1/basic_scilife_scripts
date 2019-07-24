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
from ast import literal_eval

#Keras
import tensorflow as tf
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
from tensorflow.keras.losses import mean_absolute_error, mean_squared_error, mean_absolute_percentage_error
#visualization
from tensorflow.keras.callbacks import TensorBoard
#Custom
from model_inputs import split_on_h_group, pad_cut
from lr_finder import LRFinder
from scipy import stats
import pdb
#Arguments for argparse module:
parser = argparse.ArgumentParser(description = '''A Neural Network for predicting
                                                distance between structural alignments .''')

parser.add_argument('dataframe', nargs=1, type= str,
                  default=sys.stdin, help = 'Path to np array.')

parser.add_argument('params_file', nargs=1, type= str,
                  default=sys.stdin, help = 'Path to file with net parameters')

parser.add_argument('out_dir', nargs=1, type= str,
                  default=sys.stdin, help = 'Path to output directory. Include /in end')

#from tensorflow.keras.backend import set_session
#config = tf.ConfigProto()
#config.gpu_options.allow_growth = True  # dynamically grow the memory used on the GPU
#config.log_device_placement = True  # to log device placement (on which device the operation ran)
                                    # (nothing gets printed in Jupyter, only if you run it standalone)
#sess = tf.Session(config=config)
#set_session(sess)  # set this TensorFlow session as the default session for Keras

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


def get_batch(batch_size,s="train"):
    """
    Create batch of n pairs
    """
    x1 = [*train_df['enc1']]
    x2 = [*train_df['enc2']]
    X = [x1, x2]
    scores = np.asarray(train_df["global_lddt"])

    random_numbers = np.random.choice(len(scores),size=(batch_size,),replace=False) #without replacement

    #initialize 2 empty arrays for the input image batch
    pairs=[np.zeros((batch_size, 300, 22)) for i in range(2)]

    # initialize vector for the targets
    targets=np.zeros((batch_size,))

    #Get batch data
    j = 0 #Placement in pairs and targets
    for i in random_numbers:
      r1 = np.random.choice(2) #Choose which sequence should end up in siamese 1 or 2
      r2 = np.setdiff1d([0,1], r1)[0]

      enc1 = np.eye(22)[literal_eval(X[r1][i])]
      enc2 = np.eye(22)[literal_eval(X[r2][i])]
      pairs[0][j,:,:] = pad_cut(enc1, 300, 22)
      pairs[1][j,:,:] = pad_cut(enc2, 300, 22)
      targets[j] = scores[i]
      j+=1

    return pairs, targets

def generate(batch_size, s="train"):
    """
    a generator for batches, so model.fit_generator can be used.
    """
    while True:
        pairs, targets = get_batch(batch_size,s)
        yield (pairs, targets)

#MAIN
args = parser.parse_args()
dataframe = args.dataframe[0]
params_file = args.params_file[0]
out_dir = args.out_dir[0]

#Assign data and labels
#Read df
complete_df = pd.read_csv(dataframe)
np.random.seed(2) #Set random seed - ensures same split every time
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

#Validation data
val_enc1 = [pad_cut(np.eye(22)[literal_eval(x)], 300, 22) for x in [*valid_df['enc1']]]
val_enc2 = [pad_cut(np.eye(22)[literal_eval(x)], 300, 22) for x in [*valid_df['enc2']]]
X_valid = [np.asarray(val_enc1), np.asarray(val_enc2)]
y_valid = np.asarray(valid_df['global_lddt'])
#Save validation data
np.savetxt(out_dir+'y_valid.txt', y_valid)

# #Take p first points in train
# X_train_p = [[],[]]
# y_train_p = []
# count_p = np.zeros(len(bins))
# p = 1000
# for i in range(0, len(y_train)):
#       pos = np.argmax(y_train[i])
#       if count_p[pos] <= p:
#           X_train_p[0].append(X_train[0][i])
#           X_train_p[1].append(X_train[1][i])
#           y_train_p.append(y_train[i])
#           count_p[pos]+=1

# X_train_p = [np.asarray(X_train_p[0]), np.asarray(X_train_p[1])] #convert to arrays
# y_train_p = np.asarray(y_train_p)

bins = np.expand_dims(bins, axis=0)
#Tensorboard for logging and visualization
log_name = str(time.time())
tensorboard = TensorBoard(log_dir=out_dir+log_name)


######MODEL######
#Parameters
net_params = read_net_params(params_file)
input_dim = (300,22)
num_classes = max(bins.shape)
kernel_size = 21 #google uses 21

#Variable params
num_res_blocks = int(net_params['num_res_blocks'])
base_epochs = int(net_params['base_epochs'])
finish_epochs = int(net_params['finish_epochs'])
filters = int(net_params['filters']) # Dimension of the embedding vector.
dilation_rate = int(net_params['dilation_rate'])  #dilation rate for convolutions
alpha = int(net_params['alpha'])
batch_size = 32 #int(net_params['batch_size'])
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

#Should have sum of two losses:
#1. How closely the predicted lddt matches the real one
#2. How closely the probability distribution of the bins match a gaussian distribution (kl divergence)
#Perhaps I should try to predict the deviation from the mean, since the scores are so centered around the mean.
# Add a customized layer to compute the absolute difference between the encodings
L1_layer = Lambda(lambda tensors:abs(tensors[0] - tensors[1]))
L1_distance = L1_layer([flat1, flat2])

# Add a customized layer to compute the absolute difference between the encodings
#L2_layer = Lambda(lambda tensors:keras.backend.sqrt(keras.backend.square(tensors[0] - tensors[1])))
#L2_distance = L2_layer([flat1, flat2])
#Dense final layer for classification
probabilities = Dense(num_classes, activation='softmax')(L1_distance)
bins_K = variable(value=bins)

def multiply(x):
  return tf.matmul(x, bins_K,transpose_b=True)

pred_vals = Lambda(multiply)(probabilities)
#The length of the validation data must be a multiple of batch size!
#Other wise you will have shape mismatches

#Custom loss
def bin_loss(y_true, y_pred):
  #Shold make this a log loss
        g_loss = mean_absolute_error(y_true, y_pred) #general, compare difference
	#log_g_loss = keras.backend.log(g_loss/100)
  #Gauss for loss
	#gauss = keras.backend.random_normal_variable(shape=(batch_size, 1), mean=0.7, scale=0.3) # Gaussian distribution, scale: Float, standard deviation of the normal distribution.
        kl_loss = keras.losses.kullback_leibler_divergence(y_true, y_pred) #better than comparing to gaussian
        sum_kl_loss = keras.backend.sum(kl_loss, axis =0)
        sum_g_loss = keras.backend.sum(g_loss, axis =0)
        sum_g_loss = sum_g_loss*alpha #This is basically a loss penalty


        #Normalize due to proportion
        kl_p = sum_kl_loss/(sum_g_loss+sum_kl_loss)
        g_p = sum_g_loss/(sum_g_loss+sum_kl_loss)

        sum_kl_loss = sum_kl_loss/kl_p
        sum_g_loss = sum_g_loss/g_p
        loss = sum_g_loss+sum_kl_loss
        #Scale with R? loss = loss/R - on_batch_end
  	#Normalize loss by percentage contributions: divide by contribution
  	#Write batch generator to avoid incompatibility in shapes
  	#problem at batch end due to kongruens
        return loss

#Custom validation loss
class IntervalEvaluation(Callback):
    def __init__(self, validation_data=(), interval=1):
        super(Callback, self).__init__()

        self.interval = interval
        self.X_val, self.y_val = validation_data

    def on_epoch_end(self, epoch, logs={}):
        if epoch % self.interval == 0:
            y_pred = self.model.predict(self.X_val, verbose=0)
            diff = [y_pred[i]-y_valid[i] for i in range(len(y_valid))]
            score = np.average(np.absolute(diff))
            #Pearson correlation coefficient
            R,pval = stats.pearsonr(y_valid, y_pred.flatten())
            R = np.round(R, decimals = 3)
            score = np.round(score, decimals = 3)
            print('epoch: ',epoch, ' score: ', score, ' R: ', R)

            np.savetxt(out_dir+'validpred_'+str(epoch)+'_'+str(score)+'_'+str(R)+'.txt', y_pred)

#Model: define inputs and outputs
model = Model(inputs = [in_1, in_2], outputs = pred_vals)
opt = optimizers.Adam() #remove clipnorm and add loss penalty
model.compile(loss=bin_loss,
              optimizer=opt)

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
#Checkpoint
filepath=out_dir+"weights-{epoch:02d}-.hdf5"
checkpoint = ModelCheckpoint(filepath, verbose=1, save_best_only=False, mode='max')

#Validation data
ival = IntervalEvaluation(validation_data=(X_valid, y_valid), interval=1)

#Summary of model
print(model.summary())

callbacks=[lrate, tensorboard, checkpoint, ival]
#Fit model
#Should shuffle uid1 and uid2 in X[0] vs X[1]
model.fit_generator(generate(batch_size),
            steps_per_epoch=int(len(train_df)/batch_size),
            epochs=num_epochs,
            #validation_data = [X_valid, y_valid],
            shuffle=True, #Dont feed continuously
            callbacks=callbacks)


