#! /usr/bin/env python3
# -*- coding: utf-8 -*-


import argparse
import sys
import matplotlib.pyplot as plt
import numpy as np
import glob
from sklearn.model_selection import train_test_split
from collections import Counter

import tensorflow.keras as keras
from tensorflow.keras.models import Sequential, Model
from tensorflow.keras.layers import Dense, Embedding, Flatten
from tensorflow.keras.layers import LSTM
from tensorflow.keras.layers import Reshape
from tensorflow.keras.layers import concatenate
from tensorflow.keras.backend import reshape

#import custom functions
from rnn_input import read_labels, rmsd_hot, get_encodings, get_locations, encoding_distributions, get_labels, label_distr, pad_cut
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

#Functions


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
X = [] #Save data
y = [] #Save labels


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
labels, values = zip(*counted_groups.items())
encoding_distributions('hist', values, 'Distribution of number of encoded pairs per H-group', 'Number of encoded pairs', 'log count', 10, out_dir, 'hgroups', True, [])

#Look at data distributions
encoding_distributions('hist', accessibilities, 'Distribution of Normalized surface accessibilities', '% surface accessibility', 'log count', 101, out_dir, 'acc', True, [])
encoding_distributions('bar',structures, 'Distribution of secondary structure elements', 'secondary structure', 'log count', 10, out_dir, 'str', True, ['G', 'H', 'I', 'T', 'E', 'B', 'S', 'C', '-'])
encoding_distributions('bar', letters, 'Distribution of amino acids in sequences', 'amino acid', 'log count', 22, out_dir, 'aa', True, ['A', 'R', 'N', 'D', 'C', 'E', 'Q', 'G', 'H', 'I', 'L', 'K', 'M', 'F', 'P', 'S', 'T', 'W', 'Y', 'V', 'X', '-'])
encoding_distributions('hist', seqlens, 'Distribution of sequence lengths', 'sequence length', 'count', 100, out_dir, 'seqlens', False, [])
 
#Look at info from TMalign and tree-puzzle
encoding_distributions('hist', Chains, 'Distribution of chain lengths', 'Chain length', 'count', 100, out_dir, 'chains', False, [] )
encoding_distributions('hist', Align_lens, 'Distribution of % aligned of shortest chain length' , '% aligned of shortest chain length', 'log count', 100, out_dir, 'aligned', True, [])
encoding_distributions('hist', Identities, 'Distribution of chain Identities', 'Identity (%)', 'count', 100, out_dir, 'id', False, [])
label_distr(ML_dists, rmsd_dists, 'seq_str', out_dir, 'ML AA distance', 'RMSD')

pdb.set_trace()
#Assign data and labels
X = encoding_list
y = rmsd_hot(rmsd_dists) #One-hot encode labels

#Split train data to use 80% for training and 10% for validation and 10 % for testing. 
X_train, X_valid, y_train, y_valid = train_test_split(X, y, test_size=0.2, random_state=42)
#Random state = 42 guarantees the split is the same every time. This can be both bad and good, depending on
#the selction. It makes the split consistent across changes to the network though.

X_valid, X_test, y_valid, y_test = train_test_split(X_valid, y_valid, test_size=0.5, random_state=42)



#Pad data/cut to maxlen     
X_train = pad_cut(X_train, 300)
X_valid = pad_cut(X_valid, 300)
X_test = pad_cut(X_test, 300)

pdb.set_trace()
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




######MODEL######
#Parameters
number_of_layers = 3

vocab_sizes = [23, 10, 102, 23, 10, 102] #needs to be num_unique+1 for Keras
batch_size = 10 #Number of alignments

epoch_length = int(len(X_train)/batch_size)
num_epochs = 1
forget_bias = 0.0 #Bias for LSTMs forget gate, reduce forgetting in beginning of training
num_nodes = 300
embedding_size = 10 # Dimension of the embedding vector. 1:1 ratio with hidden nodes. Should probably have smaller embeddings


init_scale = 0.1
start_learning_rate = 5.0
max_grad_norm = 0.25
keep_prob = 0.9
epsilon = 0.0000001
penalty = 1.3


    






#Define 6 different embeddings and append to model:


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
#cat_embeddings = Flatten(cat_embeddings)


lstm_out = LSTM(128, dropout=0.2, recurrent_dropout=0.2)(cat_embeddings)
outp = Dense(101, activation='sigmoid')(lstm_out)


model = Model(inputs = [embed1_in, embed2_in, embed3_in, embed4_in, embed5_in, embed6_in], outputs = outp)



#compile
model.compile(loss='binary_crossentropy',
              optimizer='adam',
              metrics=['accuracy'])


#Write summary of model
model.summary()

#Fit model
#X_train = np.transpose(X_train)



#pdb.set_trace()
model.fit(X_train, y_train,
          batch_size=batch_size,
          epochs=2,
          validation_data=(X_valid, y_valid)
          )



