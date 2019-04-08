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
from rnn_input import read_tsv, rmsd_hot
import pdb




#Arguments for argparse module:
parser = argparse.ArgumentParser(description = '''A Recurrent Neural Network for predicting
								RMSD between structural alignments based on sequences from per-residue alignments,
								secondary structure and surface accessibility.''')
 
parser.add_argument('dist_file', nargs=1, type= str,
                  default=sys.stdin, help = 'Path to distance file. Format: uid1	uid2	MLdist	RMSD')

parser.add_argument('encoding_dir', nargs=1, type= str,
                  default=sys.stdin, help = '''Path to files with encodings of alignments, secondary structure and surface acc.
                  Include /in end''')

parser.add_argument('out_dir', nargs=1, type= str,
                  default=sys.stdin, help = 'Path to output directory. Include /in end')

#Functions
#def load_data():

def plot_distr(y, name, out_dir):
	'''plot distribution of labels (normalized rmsd values)
	'''

	#Convert back to ints
	y_ints = np.argmax(y, axis = 1)

	plt.hist(y_ints, bins = 101)
	plt.savefig(out_dir+name+'.png')
	plt.close()
	
	#plt.show()
	return None


#MAIN
args = parser.parse_args()
dist_file = args.dist_file[0]
encode_dir = args.one_hot_dir[0]
out_dir = args.out_dir[0]
#Read tsv
(uids, rmsd_dists_t, rmsd_dists) = read_tsv(dist_file, 6)
#Format rmsd_dists into one-hot encoding
rmsd_dists_hot = rmsd_hot(rmsd_dists_t)


#Get macthing alignments, sendary structure and surface acc
#Should have absolute paths, will be much faster than globbing everything. Should change all_dist_rmsd.tsv to include H-family name
#raise IOerror if file is not found

X = [] #Save data
y = [] #Save labels


#Test, load only little data
max_aln_len = 0 #maximum aln length
for i in range(0,len(uids)):
	file_name = glob.glob(one_hot_dir + '*/'+uids[i]+'.enc')

	if file_name:

		matrix = np.loadtxt(file_name[0])
		if len(matrix) > max_aln_len:
			max_aln_len = len(matrix)
		X.append(matrix)
		y.append(rmsd_dists_hot[i])


#Split train data to use 80% for training and 10% for validation and 10 % for testing. 
X_train, X_valid, y_train, y_valid = train_test_split(X, y, test_size=0.2, random_state=42)
#Random state = 42 guarantees the split is the same every time. This can be both bad and good, depending on
#the selction. It makes the split consistent across changes to the network though.

X_valid, X_test, y_valid, y_test = train_test_split(X_valid, y_valid, test_size=0.5, random_state=42)

print('Train:',len(X_train), 'Valid:',len(X_valid), 'Test:',len(X_test))


#Plot distributions of labels - make sure split preserves variation in rmsd
labels = [y_train, y_valid, y_test]
names = ['y_train', 'y_valid', 'y_test']
for i in range(0,3):
	plot_distr(labels[i], names[i], out_dir)



######MODEL######
#Parameters
number_of_layers = 3
num_unrollings = max_aln_len#length of longest alignment? This should not be variable
batch_size = 2 #Number of alignments
input_size = 101
epoch_length = len(X_train)/batch_size
num_epochs = 1
forget_bias = 0.0 #Bias for LSTMs forget gate, reduce forgetting in beginning of training
num_nodes = 300
embedding_size = num_nodes # Dimension of the embedding vector. 1:1 ratio with hidden nodes
init_scale = 0.1
start_learning_rate = 5.0
max_grad_norm = 0.25
keep_prob = 0.9
num_steps = 50000
epsilon = 0.0000001
penalty = 1.3

#Define Graph
graph = tf.Graph()

with graph.as_default():
    
    
    # Classifier weights and biases. Must be the size of the number of parameters = 66 (2*33), otherwise all parameters will not be represented in logits
  softmax_w = tf.Variable(tf.random_uniform(shape = [num_nodes, input_size], minval = -init_scale, maxval = init_scale, name = 'softmax_w'))
  softmax_b = tf.Variable(tf.zeros([input_size]), name = 'softmax_b')
  
    #Embedding vector, input_size should be vocab_siz = 101
  embeddings = tf.Variable(tf.random_uniform(shape = [66, embedding_size], minval = -init_scale, maxval = init_scale), name = 'embeddings')
    


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

  logits = tf.nn.xw_plus_b(tf.concat(outputs, 0), softmax_w, softmax_b) #Computes matmul, need to have this tf concat, any other and it complains
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








pdb.set_trace()

##########RUN#########
with tf.Session(graph=graph) as session:
  tf.global_variables_initializer().run()
  print('Initialized')

  for i in range(num_epochs):
  	j = 0 #
  	for j in range(epoch_length):

  		train_feed_inputs = X_train[j_j+batch_size]
   		#Feed dict
  		feed_dict= {train_inputs: train_feed_inputs, train_labels: train_feed_labels, valid_inputs: valid_feed_inputs, valid_labels: valid_feed_labels, global_step: step, keep_probability: keep_prob}

  		_, t_perplexity, train_pred, summary = session.run([optimize, train_perplexity, train_predictions, merged], feed_dict= feed_dict)
