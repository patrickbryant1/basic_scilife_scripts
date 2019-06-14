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

import tensorflow as tf
import tensorflow.keras as keras
from tensorflow.keras import regularizers, backend
from tensorflow.keras.constraints import max_norm
from tensorflow.keras.models import Sequential, Model
from tensorflow.keras.layers import Dense, Embedding, Flatten
from tensorflow.keras.layers import Bidirectional,CuDNNLSTM, Dropout, BatchNormalization
from tensorflow.keras.layers import Reshape, Activation, RepeatVector, Permute, multiply, Lambda
from tensorflow.keras.layers import concatenate, add, Conv1D
from tensorflow.keras.callbacks import ModelCheckpoint, LearningRateScheduler, Callback
from tensorflow.keras.preprocessing.sequence import TimeseriesGenerator
from tensorflow.keras.callbacks import TensorBoard


#import custom functions
from lr_finder import LRFinder
import pdb
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

#Arguments for argparse module:
parser = argparse.ArgumentParser(description = '''A Recurrent Neural Network for predicting
                                                MSD between structural alignments .''')

parser.add_argument('dataframe', nargs=1, type= str,
                  default=sys.stdin, help = 'Path to dataframe in .csv.')

parser.add_argument('out_dir', nargs=1, type= str,
                  default=sys.stdin, help = 'Path to output directory. Include /in end')

parser.add_argument('params_file', nargs=1, type= str,
                  default=sys.stdin, help = 'Path to file with net parameters')


#MAIN
args = parser.parse_args()
dataframe = args.dataframe[0]
out_dir = args.out_dir[0]
params_file = args.params_file[0]

#Assign data and labels
#Read df
complete_df = pd.read_csv(dataframe)
#Get encodings
enc1 = []
enc2 = []
[enc1.append(literal_eval(x)) for x in complete_df['enc1']]
[enc2.append(literal_eval(x)) for x in complete_df['enc2']]
#Get longest alignment
enc_lens = []
[enc_lens.append(len(x)) for x in enc1]
#Get MLAAdist
evdist = complete_df['MLAAdist_x']
#Convert to array
X = [np.asarray(enc1), np.asarray(enc2), np.asarray(evdist)]
X = np.asarray(X)#Convert to np array
X = X.T #transpose
#One-hot encode binned data
#-4.2478745796221045 2.871030140344752
bins = np.arange(0.5,2.5,0.25)
bins = np.insert(bins,0, 0)
bins = np.append(bins, 4.5)
#Bin the TMscore RMSDs
rmsds = complete_df['RMSD_x']


binned_rmsds = np.digitize(rmsds, bins)
np.asarray(binned_rmsds)
y = np.eye(len(bins))[binned_rmsds] #deviations_hot

#Split train data to use 80% for training and 10% for validation and 10 % for testing.
X_train, X_valid, y_train, y_valid = train_test_split(X, y, test_size=0.2, random_state=42)
#Random state = 42 guarantees the split is the same every time. This can be both bad and good, depending on
#the selction. It makes the split consistent across changes to the network though.
#Get test data
X_valid, X_test, y_valid, y_test = train_test_split(X_valid, y_valid, test_size=0.5, random_state=42)


#Transpose back
X_train = X_train.T
X_valid = X_valid.T
X_test = X_test.T

#Pad data
#Get longest alignment
padlen = max(enc_lens)
#Set as two separate arrays for encoding representation
X_train_1 = pad_data(X_train[0], padlen)
X_train_2 = pad_data(X_train[1], padlen)
X_train_3 = [np.repeat(x, padlen) for x in X_train[2]] #Repeat evdist to feed at each step in lstm
X_train = [np.asarray(X_train_1),np.asarray(X_train_2)] #, np.asarray(X_train_3)]
X_valid_1 = pad_data(X_valid[0], padlen)
X_valid_2 = pad_data(X_valid[1], padlen)
X_valid_3 = [np.repeat(x, padlen) for x in X_valid[2]]
X_valid = [np.asarray(X_valid_1),np.asarray(X_valid_2)] #,  np.asarray(X_valid_3)]

#X_test = [pad_data(X_test[0], padlen),pad_data(X_test[1], padlen)]



#Tensorboard for logging and visualization
log_name = str(time.time())
tensorboard = TensorBoard(log_dir=out_dir+log_name)

######MODEL######
#Parameters
net_params = read_net_params(params_file)

#Variable params
num_res_blocks = int(net_params['num_res_blocks'])
base_epochs = int(net_params['base_epochs'])
finish_epochs = int(net_params['finish_epochs'])
num_nodes = int(net_params['num_nodes']) #Number of nodes in LSTM
embedding_size = int(net_params['embedding_size']) # Dimension of the embedding vector.
drop_rate = float(net_params['drop_rate'])  #Fraction of input units to drop
find_lr = bool(int(net_params['find_lr']))


#Fixed params
num_classes = len(bins)
vocab_sizes = [22, 22]
batch_size = 20 #Number of alignments
num_epochs = base_epochs+finish_epochs


lambda_recurrent = 0.01 #How much the model should be penalized #Lasso regression?
recurrent_max_norm = 2.0


#LR schedule
max_lr = 0.01
min_lr = max_lr/10
lr_change = (max_lr-min_lr)/(base_epochs/2-1) #Reduce further lst three epochs

#####LAYERS#####

#Define 6 different embeddings and cat
embed1_in = keras.Input(shape = [padlen])
embed1 = Embedding(vocab_sizes[0] ,embedding_size, input_length = padlen)(embed1_in)#None indicates a variable input length
#embed1 = Flatten()(embed1) #Has shape none,none,embedding_size - incompatible with evdist_in
embed2_in = keras.Input(shape =  [padlen])
embed2 = Embedding(vocab_sizes[1] ,embedding_size, input_length = padlen)(embed2_in)#None indicates a variable input length
#embed2 = Flatten()(embed2)
#evdist_in = keras.Input(shape = [padlen])

cat_embeddings = concatenate([(embed1), (embed2)]) #, (evdist_in)])
#cat_embeddings = BatchNormalization()(cat_embeddings) #Bacth normalize, focus on segment of input
#cat_embeddings = Reshape((padlen*embedding_size*2+padlen,1))(cat_embeddings)
cat_embeddings = Dropout(rate = drop_rate, name = 'cat_embed_dropped')(cat_embeddings) #Dropout


def resnet(cat_embeddings, num_res_blocks, num_classes=num_classes):
	"""Builds a resnet with bidirectinoal LSTMs of the defined depth.
	"""

	x = cat_embeddings
	# Instantiate the stack of residual units
	for res_block in range(num_res_blocks):
		lstm_out1 = Bidirectional(CuDNNLSTM(num_nodes, recurrent_regularizer = regularizers.l2(lambda_recurrent),  kernel_constraint=max_norm(recurrent_max_norm), return_sequences=True))(x) #stateful: Boolean (default False). If True, the last state for each sample at index i in a batch will be used as initial state for the sample of index i in the following batch.
		#lstm_out1 = BatchNormalization()(lstm_out1) #Bacth normalize, focus on segment of lstm_out1
		lstm_out1 = Dropout(rate = drop_rate)(lstm_out1) #Dropout
		lstm_out2 = Bidirectional(CuDNNLSTM(num_nodes, recurrent_regularizer = regularizers.l2(lambda_recurrent),  kernel_constraint=max_norm(recurrent_max_norm), return_sequences=True))(lstm_out1)
		#lstm_out2 = BatchNormalization()(lstm_out2) #Bacth normalize, focus on segment of lstm_out2
		lstm_out2 = Dropout(rate = drop_rate)(lstm_out2) #Dropout
		lstm_add1 = add([lstm_out1, lstm_out2]) #Skip connection, add before or after dropout?
		if res_block == (num_res_blocks-1): #The last block should use all time steps, not only the last output.
			lstm_out3 = Bidirectional(CuDNNLSTM(num_nodes, recurrent_regularizer = regularizers.l2(lambda_recurrent),  kernel_constraint=max_norm(recurrent_max_norm)))(lstm_add1)
		else:
			lstm_out3 = Bidirectional(CuDNNLSTM(num_nodes, recurrent_regularizer = regularizers.l2(lambda_recurrent),  kernel_constraint=max_norm(recurrent_max_norm), return_sequences=True))(lstm_add1)
		lstm_out3 = Dropout(rate = drop_rate)(lstm_out3) #Dropout
		x = add([lstm_out2, lstm_out3])

	return x


#Create resnet and get outputs
x = resnet(cat_embeddings, num_res_blocks, num_classes)



#Attention layer - information will be redistributed in the backwards pass
attention = Dense(1, activation='tanh')(x) #Normalize and extract info with tanh activated weight matrix (hidden attention weights)
attention = Flatten()(attention) #Make 1D
attention = Activation('softmax')(attention) #Softmax on all activations (normalize activations)
attention = RepeatVector(num_nodes*2)(attention) #Repeats the input "num_nodes" times.
attention = Permute([2, 1])(attention) #Permutes the dimensions of the input according to a given pattern. (permutes pos 2 and 1 of attention)

sent_representation = multiply([x, attention]) #Multiply input to attention with normalized activations
sent_representation = Lambda(lambda xin: keras.backend.sum(xin, axis=-2), output_shape=(num_nodes*2,))(sent_representation) #Sum all attentions

#Dense final layer for classification
probabilities = Dense(num_classes, activation='softmax')(sent_representation)


#Model: define inputs and outputs
model = Model(inputs = [embed1_in, embed2_in], outputs = probabilities)

#Compile model
model.compile(loss='categorical_crossentropy', #[categorical_focal_loss(alpha=.25, gamma=2)],
              optimizer='adam',
              metrics=['accuracy'])

#Write summary of model
model.summary()

#Checkpoint
filepath=out_dir+"weights-improvement-{epoch:02d}-{val_acc:.2f}.hdf5"
checkpoint = ModelCheckpoint(filepath, monitor='val_acc', verbose=1, save_best_only=False, mode='max')

#Confusion matrix
class CollectOutputAndTarget(Callback):
    def __init__(self):
        super(CollectOutputAndTarget, self).__init__()
        self.targets = []  # collect y_true batches
        self.outputs = []  # collect y_pred batches

        # the shape of these 2 variables will change according to batch shape
        # to handle the "last batch", specify `validate_shape=False`
        self.var_y_true = tf.Variable(0., validate_shape=False)
        self.var_y_pred = tf.Variable(0., validate_shape=False)

    def on_batch_end(self, batch, logs=None):
        # evaluate the variables and save them into lists
        self.targets.append(backend.eval(self.var_y_true))
        self.outputs.append(backend.eval(self.var_y_pred))
    def on_epoch_end(self, epoch, logs=None):
        y_true = []
        y_pred = []
        for i in range(len(self.targets)):
                [y_true.append(j) for j in self.targets[i]]
                [y_pred.append(j) for j in self.outputs[i]]
        y_true = np.argmax(np.array(y_true), axis = 1)
        y_pred = np.argmax(np.array(y_pred), axis = 1)
        report = classification_report(y_true, y_pred)
        with open(out_dir+'clf_report.txt', "a") as file:
                file.write('epoch: '+str(epoch)+'\n')
                file.write(report)
        self.targets = [] #reset
        self.outputs = []
# initialize the variables and the `tf.assign` ops
cbk = CollectOutputAndTarget()
fetches = [tf.assign(cbk.var_y_true, model.targets[0], validate_shape=False),
           tf.assign(cbk.var_y_pred, model.outputs[0], validate_shape=False)]
model._function_kwargs = {'fetches': fetches}  # use `model._function_kwargs` if using `Model` instead of `Sequential`


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
    # model is a Keras model
    lr_finder = LRFinder(model)

    # Train a model with batch size 20 for 5 epochs
    # with learning rate growing exponentially from 0.000001 to 1
    lr_finder.find(X_train, y_train, start_lr=0.000001, end_lr=1, batch_size=20, epochs=1)
    # Print lr and loss

    x  = lr_finder.lrs
    y = lr_finder.losses
    with open(out_dir+params_file.split('/')[-1].split('.')[0]+'.lr', "w") as file:
        for i in range(min(len(x), len(y))):
          	file.write(str(x[i]) + '\t' + str(y[i]) + '\n')
    #pdb.set_trace()
else:
  lrate = LearningRateScheduler(lr_schedule)
  callbacks=[tensorboard, checkpoint, lrate, cbk]
  validation_data=(X_valid, y_valid)





  #Fit model
  model.fit(X_train, y_train, batch_size = batch_size,
              epochs=num_epochs,
              validation_data=validation_data,
              shuffle=True, #Dont feed continuously
	      callbacks=callbacks)


  #Save model for future use

  #from tensorflow.keras.models import model_from_json
  #serialize model to JSON
  model_json = model.to_json()
  with open(out_dir+"model.json", "w") as json_file:
  	json_file.write(model_json)
