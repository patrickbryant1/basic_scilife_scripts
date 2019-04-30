#! /usr/bin/env python3
# -*- coding: utf-8 -*-


import argparse
import sys
import matplotlib.pyplot as plt
import numpy as np
import glob
from sklearn.model_selection import train_test_split
import tensorflow as tf


#import custom functions
from rnn_input import read_labels, rmsd_hot, make_dict, get_locations, get_labels, word_distributions, label_distr
import pdb




#Arguments for argparse module:
parser = argparse.ArgumentParser(description = '''A Recurrent Neural Network for predicting
								RMSD between structural alignments based on sequences from per-residue alignments,
								secondary structure and surface accessibility.''')
 
parser.add_argument('dist_file', nargs=1, type= str,
                  default=sys.stdin, help = 'Path to distance file. Format: uid1	uid2	MLdist	TMinfo..')

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


dictionary = {} #dict to save all possible combinations
accessibilities = [] #list to save all accessibilities
seqlens = [] #list to save all sequence lengths
encodings = {} #Save all encodings
#Test, load only little data
max_aln_len = 0 #maximum aln length

#Get file locations
locations = get_locations(encode_locations)

for file_name in locations:
  (encoding, dictionary, accessibilities, seqlens) = make_dict(file_name, dictionary, accessibilities, seqlens)

  file_name = file_name.split('/')
  h_group = file_name[-2]
  uid_pair = file_name[-1].split('.')[0]
  encodings[uid_pair] = encoding



pdb.set_trace()

#Get corresponding labels (rmsds) for the encoded sequences
(uids, encoding_list, rmsd_dists, ML_dists) = get_labels(encodings, distance_dict)
pdb.set_trace()

#Assign data and labels
X = np.array(encoding_list)
X = [np.array(enc) for enc in X]
y = rmsd_hot(rmsd_dists) #One-hot encode labels

#Split train data to use 80% for training and 10% for validation and 10 % for testing. 
X_train, X_valid, y_train, y_valid = train_test_split(X, y, test_size=0.2, random_state=42)
#Random state = 42 guarantees the split is the same every time. This can be both bad and good, depending on
#the selction. It makes the split consistent across changes to the network though.

X_valid, X_test, y_valid, y_test = train_test_split(X_valid, y_valid, test_size=0.5, random_state=42)

#Get all lengths for all sequences
trainlen = [len(i) for i in X_train]
validlen = [len(i) for i in X_valid]
testlen = [len(i) for i in X_test]


print('Train:',len(X_train), 'Valid:',len(X_valid), 'Test:',len(X_test))


#Plot distributions of labels and words- make sure split preserves variation in rmsd
#labels = [y_train, y_valid, y_test]
#data = (X_train, X_valid, X_test)
#names = ['Train', 'Valid', 'Test']
#for i in range(0,3):
#	label_distr(labels[i], names[i], out_dir)
#	word_distributions(data[i], 500, out_dir, names[i])



######MODEL######
#Parameters
number_of_layers = 3
num_unrollings = 0#length of longest alignment in batch
vocab_size = len(dictionary)
batch_size = 10 #Number of alignments

epoch_length = int(len(X_train)/batch_size)
num_epochs = 1
forget_bias = 0.0 #Bias for LSTMs forget gate, reduce forgetting in beginning of training
num_nodes = 300
embedding_size = num_nodes # Dimension of the embedding vector. 1:1 ratio with hidden nodes. Should probably have smaller embeddings


init_scale = 0.1
start_learning_rate = 5.0
max_grad_norm = 0.25
keep_prob = 0.9
epsilon = 0.0000001
penalty = 1.3



#Define Graph
graph = tf.Graph()

with graph.as_default():
    
    
    # Classifier weights and biases. Must be vocab_size, otherwise all words will not be represented in logits
  softmax_w = tf.Variable(tf.random_uniform(shape = [num_nodes, 101], minval = -init_scale, maxval = init_scale, name = 'softmax_w'))
  softmax_b = tf.Variable(tf.zeros([101]), name = 'softmax_b')
  
    #Embedding vector, input_size should be vocab_siz = 101
  embeddings = tf.Variable(tf.random_uniform(shape = [vocab_size, embedding_size], minval = -init_scale, maxval = init_scale), name = 'embeddings')
    

  #Number of unrollings - longest sequence length in batch
  num_unrollings = tf.placeholder(tf.uint32)
  # Input data. Create sturcutre for input data
  #Train data
  train_inputs = tf.placeholder(tf.int32, shape=[num_unrollings, batch_size])
  train_labels = tf.placeholder(tf.int32, shape=[101, batch_size])  

  #Valid data
  #valid_inputs = tf.placeholder(tf.int32, shape=[batch_size])
  #valid_labels = tf.placeholder(tf.int32, shape=[batch_size])
  #valid_hot = tf.one_hot(valid_labels, 10000)

  #Keep prob
  keep_probability = tf.placeholder(tf.float32)



  #Embed scaling
  #embed_scaling = tf.constant(1/(1-keep_prob))

  

  #Define the LSTM cell erchitecture to be used  
  def lstm_cell(keep_probability):
    cell = tf.contrib.rnn.BasicLSTMCell(num_nodes, forget_bias=forget_bias, activation = tf.nn.elu)
    return tf.contrib.rnn.DropoutWrapper(cell, output_keep_prob=keep_probability, variational_recurrent = True, dtype = tf.float32) #Only applies dropout to output weights. Can also apply to state, input, forget.
    #The variational_recurrent = True applies the same dropout mask in every step - allowing more long-term dependencies to be learned   
  
  stacked_lstm = tf.contrib.rnn.MultiRNNCell([lstm_cell(keep_probability) for _ in range(number_of_layers)])

  #LSTM
  # Initial state of the LSTM memory.
  initial_state = stacked_lstm.zero_state(batch_size, tf.float32) #Initial state. Return zero-filled state tensor(s).
  state = initial_state
  outputs = [] #Store outputs
    
  #Unrolled lstm loop  
  for i in range(num_unrollings):
    
      # The value of state is updated after processing each batch of words.
      # Look up embeddings for inputs.
      embed = tf.nn.embedding_lookup(embeddings, train_inputs[i])
     
      #Output, state of  LSTM
      output, state = stacked_lstm(embed, state)

      outputs.append(output)

  #Save final state for validation and testing
  final_state = state

  logits = tf.nn.xw_plus_b(tf.concat(outputs,0), softmax_w, softmax_b) #Computes matmul, need to have this tf concat, any other and it complains
  logits = tf.layers.batch_normalization(logits, training=True) #Batch normalize to avoid vanishing gradients
  logits = tf.reshape(logits, [num_unrollings , batch_size, input_size])   


  #Returns 1D batch-sized float Tensor: The log-perplexity for each sequence.
  #The labels are one hot encoded. 
  
  loss = tf.nn.softmax_cross_entropy_with_logits_v2(labels=tf.concat(train_labels, 0), logits=logits) 
  #Linearly constrained weights to reduce angle bias
  #true_train_perplexity = tf.math.exp(tf.reduce_mean(loss)) #Train perplexity before addition of penalty
  #loss = tf.cond(tf.abs(tf.reduce_sum(softmax_w)) > epsilon, lambda:tf.multiply(penalty, loss), lambda:tf.add(loss, 0)) #condition, TRUE, FALSE
  train_perplexity = tf.math.exp(tf.reduce_mean(loss)) #Reduce mean is a very "strong" mean
  train_predictions = tf.argmax(tf.nn.softmax(logits), axis = -1)

  
#Optimizer
  optimizer=tf.train.GradientDescentOptimizer(learning_rate=lr) #.minimize(train_perplexity) #regular SGD has been found to outpeform adaptive
  gradients, variables = zip(*optimizer.compute_gradients(train_perplexity))
  gradients, _ = tf.clip_by_global_norm(gradients, max_grad_norm)
  optimize = optimizer.apply_gradients(zip(gradients, variables))







##########RUN#########

with tf.Session(graph=graph) as session:
	tf.global_variables_initializer().run()
	print('Initialized')
	for i in range(num_epochs):

  		for j in range(epoch_length):

  			train_feed_inputs = np.array(X_train[j:j+batch_size])


  			maxlen = max(trainlen[j:j+batch_size])
  			#Pad
  			train_feed_inputs = [np.pad(inp, (0,maxlen-len(inp)), 'constant') for inp in train_feed_inputs]

  			train_feed_inputs = np.array(train_feed_inputs).T

  			train_feed_labels = y_train[j:j+batch_size]
  			pdb.set_trace()
   			#Feed dict
  			feed_dict= {train_inputs: train_feed_inputs, train_labels: train_feed_labels, keep_probability: keep_prob, num_unrollings: maxlen}

  			_, t_perplexity, train_pred, summary = session.run([optimize, train_perplexity, train_predictions, merged], feed_dict= feed_dict)


