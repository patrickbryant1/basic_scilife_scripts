#!/usr/share/python
# -*- coding: utf-8 -*-


import argparse
import sys
import numpy as np
from ast import literal_eval
import pandas as pd
import glob
from os import makedirs
from os.path import exists, join

#Preprocessing
from collections import Counter
import math
import time
from ast import literal_eval
from sklearn.model_selection import train_test_split, StratifiedShuffleSplit

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
from tensorflow.keras.losses import mean_absolute_error, mean_squared_error
#visualization
from tensorflow.keras.callbacks import TensorBoard
#Custom
from lr_finder import LRFinder
import pdb



#Arguments for argparse module:
parser = argparse.ArgumentParser(description = '''A Neural Network for predicting
                                                H-group assignment from sequence.''')

parser.add_argument('encodings', nargs=1, type= str,
                  default=sys.stdin, help = 'Path to np array with encoded aa sequences.')

parser.add_argument('df', nargs=1, type= str,
                  default=sys.stdin, help = 'Path to df with labels and other info.')

parser.add_argument('params_file', nargs=1, type= str,
                  default=sys.stdin, help = 'Path to file with net parameters')

parser.add_argument('out_dir', nargs=1, type= str,
                  default=sys.stdin, help = 'Path to output directory. Include /in end')

#Set random seed
np.random.seed(0)
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

def get_batch(batch_size,s="train"):
    """
    Create batch of n pairs, half same class, half different class
    """
    if s == 'train':
        X = X_train
        y = y_train
    else:
        X = X_test
        y = y_test

    # #n_classes, n_examples, w, h = X.shape
    #
    # # randomly sample several classes to use in the batch
    random_numbers = np.random.choice(len(y),size=(batch_size,),replace=False) #without replacement
    not_chosen = np.setdiff1d(categories,random_numbers) #not chosen categories
    # initialize 2 empty arrays for the input image batch
    pairs=[np.zeros((batch_size, 300, 21)) for i in range(2)]

    # initialize vector for the targets
    targets=np.zeros((batch_size,))
    #
    # make one half of it '1's, so 2nd half of batch has same class (=1)
    targets[batch_size//2:] = 1

    for i in range(batch_size):
        
        batch_x[i,:,:] = pad_cut(X[random_numbers[i]], 300, 21)
        batch_y[i,:] = y[random_numbers[i]]
    return batch_x, batch_y

def generate(batch_size, s="train"):
    """
    a generator for batches, so model.fit_generator can be used.
    """
    while True:
        pairs, targets = get_batch(batch_size,s)
        yield (pairs, targets)

#MAIN
args = parser.parse_args()
encodings = args.encodings[0]
df = pd.read_csv(args.df[0])
params_file = args.params_file[0]
out_dir = args.out_dir[0]

#Assign data and labels
#Read data
X = np.load(encodings, allow_pickle=True)

#Get 5 first of all above 5
groups = [*df['group_enc']]
counted_groups = Counter(groups)
max5labels = []
max5onehot = []
for group in [*counted_groups.keys()]:
    if counted_groups[group]>5:
        ind = np.asarray(df[df['group_enc'] == group].index)
        selected = np.random.choice(ind, 5, replace = False)
        for i in selected:
            max5labels.append(group)
            max5onehot.append(X[i])
    else:
        ind = np.asarray(df[df['group_enc'] == group].index)
        for i in ind:
            max5labels.append(group)
            max5onehot.append(X[i])

#Create np arrays
X = np.asarray(max5onehot)
y = np.asarray(max5labels)

#Split so both groups are represented in train and test
sss = StratifiedShuffleSplit(n_splits=1, test_size=0.4, random_state=0)
for train_index, test_index in sss.split(X, y):
    X_train, X_valid = X[train_index], X[test_index]
    y_train, y_valid = y[train_index], y[test_index]

#Onehot encode labels
num_classes = max(y)+1
y_train = np.eye(max(y)+1)[y_train]
y_valid = np.eye(max(y)+1)[y_valid]

#Pad X_valid
padded_X_valid = []
for i in range(0,len(X_valid)):
    padded_X_valid.append(pad_cut(X_valid[i], 300, 21))
X_valid = np.asarray(padded_X_valid)

######MODEL######
#Parameters
net_params = read_net_params(params_file)
batch_size = 16
input_dim = (300,21)
num_classes = num_classes
kernel_size = 21 #google uses 21
seq_length=300
#Variable params
num_res_blocks = 1#int(net_params['num_res_blocks'])
base_epochs = int(net_params['base_epochs'])
finish_epochs = int(net_params['finish_epochs'])
filters = int(net_params['filters']) # Dimension of the embedding vector.
dilation_rate = int(net_params['dilation_rate'])  #dilation rate for convolutions

#lr opt
find_lr = False
#LR schedule
step_size = 5
num_cycles = 3
num_epochs = step_size*2*num_cycles
num_steps = int(len(y_train)/batch_size)
max_lr = 0.0009
min_lr = max_lr/10
lr_change = (max_lr-min_lr)/step_size  #(step_size*num_steps) #How mauch to change each batch
lrate = min_lr

#MODEL
in_1 = keras.Input(shape = input_dim)
#in_2 = keras.Input(shape = input_dim)
def resnet(x, num_res_blocks):
	"""Builds a resnet with 1D convolutions of the defined depth.
	"""


    	# Instantiate the stack of residual units
    	#Similar to ProtCNN, but they used batch_size = 64, 2000 filters and kernel size of 21
	for res_block in range(num_res_blocks):
		batch_out1 = BatchNormalization()(x) #Bacth normalize, focus on segment
		activation1 = Activation('relu')(x)
		conv_out1 = Conv1D(filters = filters, kernel_size = kernel_size, dilation_rate = dilation_rate, input_shape=input_dim, padding ="same")(activation1)
		batch_out2 = BatchNormalization()(conv_out1) #Bacth normalize, focus on segment
		activation2 = Activation('relu')(conv_out1)
        #Downsample - half filters
		conv_out2 = Conv1D(filters = int(filters/2), kernel_size = kernel_size, dilation_rate = 1, input_shape=input_dim, padding ="same")(activation2)
		x = Conv1D(filters = int(filters/2), kernel_size = kernel_size, dilation_rate = 1, input_shape=input_dim, padding ="same")(x)
		x = add([x, conv_out2]) #Skip connection



	return x

#Initial convolution
in_1_conv = Conv1D(filters = filters, kernel_size = kernel_size, dilation_rate = 2, input_shape=input_dim, padding ="same")(in_1)
#in_2_conv = Conv1D(filters = filters, kernel_size = kernel_size, dilation_rate = 2, input_shape=input_dim, padding ="same")(in_2)
#Output (batch, steps(len), filters), filters = channels in next
x1 = resnet(in_1_conv, num_res_blocks)
#x2 = resnet(in_2_conv, num_res_blocks)
drop1 = Dropout(0.8)(x1) #Fraction to drop. In Keras dropout is disabled in test mode
#Average pool along sequence axis
#x = AveragePooling1D(data_format='channels_first')(x) #data_format='channels_first'
maxpool1 = MaxPooling1D(pool_size=seq_length)(drop1)
#maxpool2 = MaxPooling1D(pool_size=seq_length)(x2)
#cat = concatenate([maxpool1, maxpool2]) #Cat convolutions

flat1 = Flatten(name='features')(maxpool1)  #Flatten
#flat2 = Flatten()(maxpool2)  #Flatten

# Add a customized layer to compute the absolute difference between the encodings
#L1_layer = Lambda(lambda tensors:abs(tensors[0] - tensors[1]))
#L1_distance = L1_layer([flat1, flat2])

# Add a customized layer to compute the absolute difference between the encodings
#L2_layer = Lambda(lambda tensors:keras.backend.sqrt(keras.backend.square(tensors[0] - tensors[1])))
#L2_distance = L2_layer([flat1, flat2])
#Dense final layer for classification
probability = Dense(num_classes, activation='sigmoid')(flat1)


#Model: define inputs and outputs
model = Model(inputs = [in_1], outputs = probability)
opt = optimizers.Adam(clipnorm=1.)
model.compile(loss='categorical_crossentropy',
              metrics = ['accuracy'],
              optimizer=opt)

#LRFinder
if find_lr == True:
  lr_finder = LRFinder(model)
  #data
  padded_X_train = []
  for i in range(0,len(X_train)):
      padded_X_train.append(pad_cut(X_train[i], 300, 21))

  lr_finder.find(np.asarray(padded_X_train), y_train, start_lr=0.001, end_lr=1, batch_size=batch_size, epochs=1)
  losses = lr_finder.losses
  lrs = lr_finder.lrs
  l_l = np.asarray([lrs, losses])
  np.savetxt(out_dir+'lrs_losses.txt', l_l)
  num_epochs = 0

#LR schedule
class LRschedule(Callback):
  '''lr scheduel according to one-cycle policy.
  '''
  def __init__(self, interval=1):
    super(Callback, self).__init__()
    self.lr_change = lr_change #How mauch to change each batch
    self.lr = min_lr
    self.interval = interval

  def on_epoch_end(self, epoch, logs={}):
    if epoch > 0 and epoch%step_size == 0:
      self.lr_change = self.lr_change*-1 #Change decrease/increase
    self.lr = self.lr + self.lr_change
    layer_name = 'features'
    intermediate_layer_model = Model(inputs=model.input,
                                 outputs=model.get_layer(layer_name).output)
    intermediate_output = np.asarray(intermediate_layer_model.predict(X_valid))
    np.save(out_dir+'emb_'+str(epoch)+'.npy', intermediate_output)

    #Set lr
    print(' ',self.lr)
    keras.backend.set_value(self.model.optimizer.lr, self.lr)

#Save y_valid
np.save(out_dir+'y_valid.npy', np.argmax(y_valid, axis = 1))
#Lrate
lrate = LRschedule()

#Tensorboard for logging and visualization
# save class labels to disk to color data points in TensorBoard accordingly
log_name = str(time.time())
log_dir=out_dir+log_name

tensorboard = TensorBoard(log_dir=log_dir)


#Checkpoint
filepath=out_dir+"weights-{epoch:02d}.hdf5"
checkpoint = ModelCheckpoint(filepath, monitor='val_acc', verbose=1, save_best_only=False, mode='max')

#Summary of model
print(model.summary())

callbacks=[lrate, checkpoint, tensorboard]
#Fit model
#Should shuffle uid1 and uid2 in X[0] vs X[1]
model.fit_generator(generate(batch_size),
             steps_per_epoch=num_steps,
             epochs=num_epochs,
             validation_data = [X_valid, y_valid],
             shuffle=False, #Dont feed continuously
             callbacks=callbacks)

 #from tensorflow.keras.models import model_from_json
 #serialize model to JSON
model_json = model.to_json()
with open(out_dir+"model.json", "w") as json_file:
     json_file.write(model_json)
