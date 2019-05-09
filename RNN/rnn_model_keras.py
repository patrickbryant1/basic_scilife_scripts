#! /usr/bin/env python3
# -*- coding: utf-8 -*-


import argparse
import sys
#import matplotlib.pyplot as plt
import numpy as np
import glob
from sklearn.model_selection import train_test_split
from sklearn.metrics import classification_report, confusion_matrix
from collections import Counter
import math
import time


import tensorflow.keras as keras
from tensorflow.keras import regularizers
from tensorflow.keras.constraints import max_norm
from tensorflow.keras.models import Sequential, Model
from tensorflow.keras.layers import Dense, Embedding, Flatten
from tensorflow.keras.layers import Bidirectional,CuDNNLSTM, Dropout
from tensorflow.keras.layers import Reshape
from tensorflow.keras.layers import concatenate
from tensorflow.keras.callbacks import ModelCheckpoint, LearningRateScheduler
from tensorflow.keras.preprocessing.sequence import TimeseriesGenerator
from tensorflow.keras.callbacks import TensorBoard

#import custom functions
from rnn_input import read_labels, rmsd_hot, get_encodings, get_locations, encoding_distributions, get_labels, label_distr, split_on_h_group, pad_cut
from lr_finder import LRFinder
import pdb




#Arguments for argparse module:
parser = argparse.ArgumentParser(description = '''A Recurrent Neural Network for predicting
                                                                RMSD between structural alignments based on sequences from per-residue alignments,
                                                                secondary structure and surface accessibility.''')
 
parser.add_argument('dist_file', nargs=1, type= str,
                  default=sys.stdin, help = 'Path to distance file. Format: uid1        uid2    MLdist  TMinfo..')

parser.add_argument('encode_locations', nargs=1, type= str,
                  default=sys.stdin, help = '''Paths to files with encodings of alignments, secondary structure and surface acc.
                  Include /in end''')

parser.add_argument('out_dir', nargs=1, type= str,
                  default=sys.stdin, help = 'Path to output directory. Include /in end')


#MAIN
args = parser.parse_args()
dist_file = args.dist_file[0]
encode_locations = args.encode_locations[0]
out_dir = args.out_dir[0]
#Read tsv
(distance_dict) = read_labels(dist_file)
#Format rmsd_dists into one-hot encoding
#rmsd_dists_hot = rmsd_hot(rmsd_dists_t)


#Get macthing alignments, sendary structure and surface acc
accessibilities = [] #list to save all accessibilities
structures = [] #list to save all secondary structures
letters = [] #List to save all amino acids
seqlens = [] #list to save all sequence lengths

encodings = {} #Save all encodings
H_groups = {} #Save all H-groups for split
threshold = 6 #ML seq distance threshold
#Test, load only little data
max_aln_len = 0 #maximum aln length

#Get file locations
locations = get_locations(encode_locations)

for file_name in locations:
  (encoding) = get_encodings(file_name)

  file_name = file_name.split('/')
  H_group = file_name[-2]
  uid_pair = file_name[-1].split('.')[0]
  encodings[uid_pair] = encoding
  H_groups[uid_pair] = H_group


#Get corresponding labels (rmsds) for the encoded sequences
(uids, encoding_list, rmsd_dists, ML_dists, Chains, Align_lens, Identities, letters, structures, accessibilities, seqlens, H_group_list) = get_labels(encodings, distance_dict, threshold, H_groups)

#Look at H-groups
counted_groups = Counter(H_group_list)
unique_groups, counts = zip(*counted_groups.items())
#encoding_distributions('hist', values, 'Distribution of number of encoded pairs per H-group', 'Number of encoded pairs', 'log count', 10, out_dir, 'hgroups', True, [])

#Look at data distributions
#encoding_distributions('hist', accessibilities, 'Distribution of Normalized surface accessibilities', '% surface accessibility', 'log count', 101, out_dir, 'acc', True, [])
#encoding_distributions('bar',structures, 'Distribution of secondary structure elements', 'secondary structure', 'log count', 10, out_dir, 'str', True, ['G', 'H', 'I', 'T', 'E', 'B', 'S', 'C', '-'])
#encoding_distributions('bar', letters, 'Distribution of amino acids in sequences', 'amino acid', 'log count', 22, out_dir, 'aa', True, ['A', 'R', 'N', 'D', 'C', 'E', 'Q', 'G', 'H', 'I', 'L', 'K', 'M', 'F', 'P', 'S', 'T', 'W', 'Y', 'V', 'X', '-'])
#encoding_distributions('hist', seqlens, 'Distribution of sequence lengths', 'sequence length', 'count', 100, out_dir, 'seqlens', False, [])
 
#Look at info from TMalign and tree-puzzle
#encoding_distributions('hist', Chains, 'Distribution of chain lengths', 'Chain length', 'count', 100, out_dir, 'chains', False, [] )
#encoding_distributions('hist', Align_lens, 'Distribution of % aligned of shortest chain length' , '% aligned of shortest chain length', 'log count', 100, out_dir, 'aligned', True, [])
#encoding_distributions('hist', Identities, 'Distribution of chain Identities', 'Identity (%)', 'count', 100, out_dir, 'id', False, [])
#label_distr(ML_dists, rmsd_dists, 'seq_str', out_dir, 'ML AA distance', 'RMSD')


#Assign data and labels
y = rmsd_hot(rmsd_dists, [0,20,40,60,80,100]) #One-hot encode labels


(X_train, y_train, X_valid, y_valid, X_test, y_test) = split_on_h_group(encoding_list, H_group_list, unique_groups, counted_groups, [0.8, 0.1, 0.1], y, out_dir)

#Split train data to use 80% for training and 10% for validation and 10 % for testing. 
#X_train, X_valid, y_train, y_valid = train_test_split(X, y, test_size=0.2, random_state=42)
#Random state = 42 guarantees the split is the same every time. This can be both bad and good, depending on
#the selction. It makes the split consistent across changes to the network though.

#X_valid, X_test, y_valid, y_test = train_test_split(X_valid, y_valid, test_size=0.5, random_state=42)



#Pad data/cut to maxlen     
X_train = pad_cut(X_train, 300)
X_valid = pad_cut(X_valid, 300)
X_test = pad_cut(X_test, 300)

#Get all lengths for all sequences
#trainlen = [len(i[0]) for i in X_train]
#validlen = [len(i[0]) for i in X_valid]
#testlen = [len(i[0]) for i in X_test]


print('Train:',len(X_train[0]), 'Valid:',len(X_valid[0]), 'Test:',len(X_test[0]))


#Plot distributions of labels and words- make sure split preserves variation in rmsd
#labels = [y_train, y_valid, y_test]
#data = (X_train, X_valid, X_test)
#names = ['Train', 'Valid', 'Test']
#for i in range(0,3):
#       pdb.set_trace()
#       (labels[i], names[i], out_dir)


#Tensorboard for logging and visualization
log_name = str(time.time())
tensorboard = TensorBoard(log_dir=out_dir+log_name)

######MODEL######
#Parameters
num_classes = 6

vocab_sizes = [23, 10, 102, 23, 10, 102] #needs to be num_unique+1 for Keras
batch_size = 20 #Number of alignments

epoch_length = int(len(X_train)/batch_size)
num_epochs = 23
forget_bias = 0.0 #Bias for LSTMs forget gate, reduce forgetting in beginning of training
num_nodes = 128 #Number of nodes in LSTM
embedding_size = 10 # Dimension of the embedding vector.
drop_rate = 0.5

lambda_recurrent = 0.01 #High much the model should be penalized #Lasso regression?
recurrent_max_norm = 2.0
    

#Opt
find_lr = False
max_lr = 0.01
min_lr = max_lr/10
lr_change = (max_lr-min_lr)/5 #Reduce further lst three epochs


#####LAYERS#####

#Define 6 different embeddings and cat
embed1_in = keras.Input(shape = [None])
embed1 = Embedding(vocab_sizes[0] ,embedding_size, input_length = None)(embed1_in)#None indicates a variable input length
embed2_in =  keras.Input(shape =  [None])
embed2 = Embedding(vocab_sizes[1] ,embedding_size, input_length = None)(embed2_in)#None indicates a variable input length
embed3_in =  keras.Input(shape =  [None])
embed3 = Embedding(vocab_sizes[2] ,embedding_size, input_length = None)(embed3_in)#None indicates a variable input length
embed4_in =  keras.Input(shape =  [None])
embed4 = Embedding(vocab_sizes[3] ,embedding_size, input_length = None)(embed4_in)#None indicates a variable input length
embed5_in =  keras.Input(shape =  [None])
embed5 = Embedding(vocab_sizes[4] ,embedding_size, input_length = None)(embed5_in)#None indicates a variable input length
embed6_in =  keras.Input(shape =  [None])
embed6 = Embedding(vocab_sizes[5] ,embedding_size, input_length = None)(embed6_in)#None indicates a variable input length

cat_embeddings = concatenate([(embed1), (embed2), (embed3), (embed4), (embed5), (embed6)])

cat_embeddings = Dropout(rate = drop_rate)(cat_embeddings) #Dropout
lstm_out1 = Bidirectional(CuDNNLSTM(num_nodes, recurrent_regularizer = regularizers.l2(lambda_recurrent),  kernel_constraint=max_norm(recurrent_max_norm), return_sequences=True))(cat_embeddings) #stateful: Boolean (default False). If True, the last state for each sample at index i in a batch will be used as initial state for the sample of index i in the following batch.
lstm_out1 = Dropout(rate = drop_rate)(lstm_out1) #Dropout
#lstm_out1 = Reshape(-1, num_nodes*2, 1)(lstm_out1) #num_nodes*2 since bidirectional LSTM
lstm_out2 = Bidirectional(CuDNNLSTM(int(num_nodes/2)))(lstm_out1)
lstm_out2 = Dropout(rate = drop_rate)(lstm_out2) #Dropout

outp = Dense(num_classes, activation='relu')(lstm_out2)


model = Model(inputs = [embed1_in, embed2_in, embed3_in, embed4_in, embed5_in, embed6_in], outputs = outp)

#compile
model.compile(loss='categorical_crossentropy',
              optimizer='adam',
              metrics=['accuracy'])


#Write summary of model
model.summary()

#Checkpoint
filepath="weights-improvement-{epoch:02d}-{val_acc:.2f}.hdf5"
checkpoint = ModelCheckpoint(filepath, monitor='val_acc', verbose=1, save_best_only=True, mode='max')


#LR schedule
def lr_schedule(epochs):
  '''lr scheduel according to one-cycle policy.
  '''
  
  #Increase lrate in beginning
  if epochs == 0:
    lrate = min_lr
  elif (epochs <6 and epochs > 0):
    lrate = min_lr+(epochs*lr_change)
  #Decrease further below min_lr last three epochs
  elif epochs > 10:
    lrate = min_lr/(10*(epochs-10))
  #After the max lrate is reached, decrease it back to the min
  else:
    lrate = max_lr-((epochs-5)*lr_change)

  print(epochs,lrate)
  return lrate


if find_lr == True:
    # model is a Keras model
    lr_finder = LRFinder(model)

    # Train a model with batch size 20 for 5 epochs
    # with learning rate growing exponentially from 0.000001 to 1
    lr_finder.find(X_train, y_train, start_lr=0.000001, end_lr=1, batch_size=20, epochs=5)
    # Plot the loss, ignore 20 batches in the beginning and 5 in the end
    lr_finder.plot_loss(n_skip_beginning=20, n_skip_end=5)

else:
  lrate = LearningRateScheduler(lr_schedule)
  callbacks=[tensorboard, checkpoint, lrate]
  steps_per_epoch = (len(X_train) / batch_size)
  validation_data=(X_valid, y_valid)


#Fit model
model.fit(X_train, y_train, batch_size = batch_size,             
              epochs=num_epochs,
              validation_data=validation_data,
              shuffle=True, #Dont feed continuously
	      callbacks=callbacks) #, lr_scheduler])


pdb.set_trace()

#Save model for future use

#from tensorflow.keras.models import model_from_json   
#serialize model to JSON
model_json = model.to_json()
with open(out_dir+"model.json", "w") as json_file:
	json_file.write(model_json)

# serialize weights to HDF5
model.save_weights(out_dir+"model.h5")
print("Saved model to disk")
