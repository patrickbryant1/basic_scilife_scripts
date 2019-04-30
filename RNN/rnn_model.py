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
from rnn_input import read_labels, rmsd_hot, get_encodings, get_locations, encoding_distributions, get_labels, label_distr
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


accessibilities = [] #list to save all accessibilities
structures = [] #list to save all secondary structures
letters = [] #List to save all amino acids
seqlens = [] #list to save all sequence lengths

encodings = {} #Save all encodings
threshold = 6 #ML seq distance threshold
#Test, load only little data
max_aln_len = 0 #maximum aln length

#Get file locations
locations = get_locations(encode_locations)

for file_name in locations:
  (encoding) = get_encodings(file_name)

  file_name = file_name.split('/')
  h_group = file_name[-2]
  uid_pair = file_name[-1].split('.')[0]
  encodings[uid_pair] = encoding


#Get corresponding labels (rmsds) for the encoded sequences
(uids, encoding_list, rmsd_dists, ML_dists, Chains, Align_lens, Identities, letters, structures, accessibilities, seqlens) = get_labels(encodings, distance_dict, threshold)

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
X = np.array(encoding_list)
y = rmsd_hot(rmsd_dists) #One-hot encode labels

#Split train data to use 80% for training and 10% for validation and 10 % for testing. 
X_train, X_valid, y_train, y_valid = train_test_split(X, y, test_size=0.2, random_state=42)
#Random state = 42 guarantees the split is the same every time. This can be both bad and good, depending on
#the selction. It makes the split consistent across changes to the network though.

X_valid, X_test, y_valid, y_test = train_test_split(X_valid, y_valid, test_size=0.5, random_state=42)

#Get all lengths for all sequences
trainlen = [len(i[0]) for i in X_train]
validlen = [len(i[0]) for i in X_valid]
testlen = [len(i[0]) for i in X_test]


print('Train:',len(X_train), 'Valid:',len(X_valid), 'Test:',len(X_test))


#Plot distributions of labels and words- make sure split preserves variation in rmsd
#labels = [y_train, y_valid, y_test]
#data = (X_train, X_valid, X_test)
#names = ['Train', 'Valid', 'Test']
#for i in range(0,3):
#	pdb.set_trace()
#	(labels[i], names[i], out_dir)




######MODEL######
#Parameters
number_of_layers = 3
num_unrollings = 0#length of longest alignment in batch
vocab_sizes = [22, 9, 101, 22, 9, 101]
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



#Define Graph
graph = tf.Graph()

with graph.as_default():
    
    
    # Classifier weights and biases. Must be vocab_size, otherwise all words will not be represented in logits
  softmax_w = tf.Variable(tf.random_uniform(shape = [num_nodes, 101], minval = -init_scale, maxval = init_scale, name = 'softmax_w'))
  softmax_b = tf.Variable(tf.zeros([101]), name = 'softmax_b')
  
    #Embedding vectors, input_size should be vocab_size (variable btq 22, 9 and 101)
    #One embedding vector for each aa, 2ndarystr, acc in each pair 
  embedding1 = tf.Variable(tf.random_uniform(shape = [vocab_sizes[0], embedding_size], minval = -init_scale, maxval = init_scale), name = 'embedding1')
  embedding2 = tf.Variable(tf.random_uniform(shape = [vocab_sizes[1], embedding_size], minval = -init_scale, maxval = init_scale), name = 'embedding2')
  embedding3 = tf.Variable(tf.random_uniform(shape = [vocab_sizes[2], embedding_size], minval = -init_scale, maxval = init_scale), name = 'embedding3')
  embedding4 = tf.Variable(tf.random_uniform(shape = [vocab_sizes[3], embedding_size], minval = -init_scale, maxval = init_scale), name = 'embedding4')
  embedding5 = tf.Variable(tf.random_uniform(shape = [vocab_sizes[4], embedding_size], minval = -init_scale, maxval = init_scale), name = 'embedding5')
  embedding6 = tf.Variable(tf.random_uniform(shape = [vocab_sizes[5], embedding_size], minval = -init_scale, maxval = init_scale), name = 'embedding6')
    

  #Number of unrollings - longest sequence length in batch
  num_unrollings = tf.placeholder(dtype = tf.uint32, shape = (1))
  # Input data. Create sturcutre for input data
  #Train data
  train_inputs = tf.placeholder(tf.int32, shape=[None, batch_size]) #None for varying input size = number of unrollings
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
  #Have to figure out how to loop over the tensor
  for i in range(>tensor length<):#range(num_unrollings): #Should go through sequence one position (= "word") at a time
    
      # The value of state is updated after processing each batch of words.
      # Look up embeddings for inputs.
      embed1 = tf.nn.embedding_lookup(embedding1, train_inputs[i][0]) #input: 1xlength of longest sequence in batch 
      embed2 = tf.nn.embedding_lookup(embedding2, train_inputs[i][1])
      embed3 = tf.nn.embedding_lookup(embedding3, train_inputs[i][2])
      embed4 = tf.nn.embedding_lookup(embedding1, train_inputs[i][3])
      embed5 = tf.nn.embedding_lookup(embedding2, train_inputs[i][4])
      embed6 = tf.nn.embedding_lookup(embedding3, train_inputs[i][5])
     
      #Reshape to one dimension
      flat1 = tf.reshape(embed1, [-1])
      flat2 = tf.reshape(embed2, [-1])
      flat3 = tf.reshape(embed3, [-1])
      flat4 = tf.reshape(embed4, [-1])
      flat5 = tf.reshape(embed5, [-1])
      flat6 = tf.reshape(embed6, [-1])

      #Cat all embeddings
      cat_embed = tf.concat([flat1, flat2, flat3, flat4, flat5, flat6],0)

      #Output, state of  LSTM
      output, state = stacked_lstm(cat_embed, state)

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


