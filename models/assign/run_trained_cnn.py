#!/usr/share/python3
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
from tensorflow.keras.models import model_from_json
import tensorflow as tf
from tensorflow.keras import regularizers,optimizers
import tensorflow.keras as keras
from tensorflow.keras.constraints import max_norm
from tensorflow.keras.callbacks import ModelCheckpoint, LearningRateScheduler, Callback
from tensorflow.keras.models import Sequential, Model
from tensorflow.keras.layers import Dense, Dropout, Activation, Conv1D, Reshape, MaxPooling1D, Dot, Masking
from tensorflow.keras.layers import Activation, RepeatVector, Permute, Multiply, Lambda, GlobalAveragePooling1D
from tensorflow.keras.layers import concatenate, add, Conv1D, BatchNormalization, Flatten, Subtract
from tensorflow.keras.backend import epsilon, clip, get_value, set_value, transpose, variable, square
from tensorflow.layers import AveragePooling1D
from tensorflow.keras.losses import mean_absolute_error, mean_squared_error

import pdb


#Arguments for argparse module:
parser = argparse.ArgumentParser(description = '''A program that reads a keras model from a .json and a .h5 file''')

parser.add_argument('json_file', nargs=1, type= str,
                  default=sys.stdin, help = 'path to .json file with keras model to be opened')

parser.add_argument('weights', nargs=1, type= str,
                  default=sys.stdin, help = '''path to .h5 file containing weights for net.''')

parser.add_argument('encodings_train', nargs=1, type= str,
                  default=sys.stdin, help = 'Path to np array with encoded aa sequences used to train the net.')

parser.add_argument('dataframe_train', nargs=1, type= str,
                  default=sys.stdin, help = '''path to dataframe with labels used to train the net.''')

parser.add_argument('encodings_test', nargs=1, type= str,
                  default=sys.stdin, help = 'Path to np array with encoded aa sequences to be used for testing.')

parser.add_argument('dataframe_test', nargs=1, type= str,
                  default=sys.stdin, help = '''path to dataframe with labels to be used for testing.''')

parser.add_argument('out_dir', nargs=1, type= str,
                  default=sys.stdin, help = '''path to output directory.''')


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

def load_model(json_file, weights):

	global model

	json_file = open(json_file, 'r')
	model_json = json_file.read()
	model = model_from_json(model_json)
	model.load_weights(weights)
	model._make_predict_function()
	#model.compile(optimizer='adam', loss='categorical_crossentropy', metrics=['accuracy'])
	return model

#MAIN
args = parser.parse_args()
json_file = (args.json_file[0])
weights = (args.weights[0])
encodings_train = args.encodings_train[0]
dataframe_train = args.dataframe_train[0]
encodings_test = args.encodings_test[0]
dataframe_test = args.dataframe_test[0]
out_dir = args.out_dir[0]

#Assign data and labels
#Read df
df = pd.read_csv(dataframe)

#Assign data and labels
#Read data
X = np.load(encodings, allow_pickle=True)
y = np.asarray([*df['group_enc']])

#Pad X
padded_X = []
for i in range(0,len(X)):
    padded_X.append(pad_cut(X[i], 300, 21))
X = np.asarray(padded_X)


#Load and run model
model = load_model(json_file, weights)

#Get embedding layers
features1 = Model(inputs=model.input, outputs=model.get_layer('features1').output)
features2 = Model(inputs=model.input, outputs=model.get_layer('features2').output)

#Get average embeddings for all entries
emb1 = np.asarray(features1.predict([X,X]))
emb2 = np.asarray(features2.predict([X,X]))
average_emb = np.average([emb1, emb2], axis = 0)
#Save embeddings
#np.save(out_dir+'average_feature_emb.npy', average_emb)

#Compute class labels by averaging the embeddings for each H-group.
class_embeddings = []
for group in range(0,max(y)+1):
    emb_match = average_emb[np.where(y == group)]
    class_emb = np.average(emb_match, axis = 0)
    class_embeddings.append(class_emb)

class_embeddings = np.asarray(class_embeddings)
#Save class embeddings
np.save(out_dir+'class_emb.npy', class_embeddings)



#Test ZSL against class embeddings
#Set random seed - use the same as wehen training the siamese net to be able to evaluate the GZSL vs SZL
np.random.seed(0)
#Split without overlaps of classes between train and test
train_classes = np.random.choice(max(y), size = (int((max(y)+1)*0.8),), replace = False)
valid_classes = np.setdiff1d(range(max(y)),train_classes) #categories not chosen for training

# #Get 5 first of all above 5
# groups = [*df['group_enc']]
# counted_groups = Counter(groups)
# max5labels = []
# max5emb1 = []
# max5emb2 = []
# max5_avemb = []
# for group in [*counted_groups.keys()]:
#     if counted_groups[group]>5:
#         ind = np.asarray(df[df['group_enc'] == group].index)
#         selected = np.random.choice(ind, 5, replace = False)
#         for i in selected:
#             max5labels.append(group)
#             max5emb1.append(emb1[i])
#             max5emb2.append(emb2[i])
#             max5_avemb.append(average_emb[i])
#
#     else:
#         ind = np.asarray(df[df['group_enc'] == group].index)
#         for i in ind:
#             max5labels.append(group)
#             max5emb1.append(emb1[i])
#             max5emb2.append(emb2[i])
#             max5_avemb.append(average_emb[i])

# emb1 = np.asarray(max5emb1)
# emb2 = np.asarray(max5emb2)
# average_emb = np.asarray(max5_avemb)
# y = np.asarray(max5labels)
# np.save(out_dir+'max5labels.npy', y)
train_index = np.isin(y, train_classes)
valid_index = np.isin(y, valid_classes)
pdb.set_trace()

def zsl_test(indices, type, out_dir):
    "A function that runs ZSL according to provided data"
    trials= [emb1[indices], emb2[indices], average_emb[indices]]
    targets = y[indices]
    names = ['emb1', 'emb2', 'average_emb']
    for t in range(3):
        item = trials[t]
        #Compute L1 distance to all class_embeddings
        pred_ranks = []
        for i in range(0, len(item)):
            true = targets[i]
            diff = np.absolute(class_embeddings-item[i])
            dists = np.sum(diff, axis = 1)
            ranks = np.argsort(dists)
            rank = np.where(ranks == true)[0][0]
            pred_ranks.append(rank)

        #Save predicted_ranks
        pred_ranks = np.asarray(pred_ranks)
        np.save(out_dir+type+'_'+names[t]+'_pred_ranks.npy', pred_ranks)

    return None

zsl_test(train_index, 'train', out_dir)
zsl_test(valid_index, 'valid', out_dir)
pdb.set_trace()
