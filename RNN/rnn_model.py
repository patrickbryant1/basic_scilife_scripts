#! /usr/bin/env python3
# -*- coding: utf-8 -*-


import argparse
import sys
import matplotlib.pyplot as plt
import numpy as np
import glob
import pandas as pd
from sklearn.model_selection import train_test_split

#import custom functions
from rnn_input import read_tsv, rmsd_hot
import pdb




#Arguments for argparse module:
parser = argparse.ArgumentParser(description = '''A Recurrent Neural Network for predicting
								RMSD between structural alignments based on sequences from per-residue alignments,
								secondary structure and surface accessibility.''')
 
parser.add_argument('dist_file', nargs=1, type= str,
                  default=sys.stdin, help = 'Path to distance file. Format: uid1	uid2	MLdist	RMSD')

parser.add_argument('one_hot_dir', nargs=1, type= str,
                  default=sys.stdin, help = 'Path to files with one-hot encodings of alignments, sendary structure and surface acc.')

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









######MODEL######


#MAIN
args = parser.parse_args()
dist_file = args.dist_file[0]
one_hot_dir = args.one_hot_dir[0]
out_dir = args.out_dir[0]
#Read tsv
(uids, rmsd_dists) = read_tsv(dist_file, 6)
#Format rmsd_dists into one-hot encoding
rmsd_dists_hot = rmsd_hot(rmsd_dists)


#Get macthing alignments, sendary structure and surface acc
#Should have absolute paths, will be much faster than globbing everything. Should change all_dist_rmsd.tsv to include H-family name
#raise IOerror if file is not found

X = [] #Save data
y = [] #Save labels


#Test, load only little data
for i in range(0,len(uids)):
	file_name = glob.glob(one_hot_dir + '*/'+uids[i]+'.hot')

	if file_name:

		matrix = np.loadtxt(file_name[0])
		X.append(matrix)
		y.append(rmsd_dists_hot[i])


#Split train data to use 80% for training and 10% for validation and 10 % for testing. 
X_train, X_valid, y_train, y_valid = train_test_split(X, y, test_size=0.2, random_state=42)
#Random state = 42 guarantees the split is the same every time. This can be both bad and good, depending on
#the selction. It makes the split consistent across changes to the network though.

X_valid, X_test, y_valid, y_test = train_test_split(X_valid, y_valid, test_size=0.5, random_state=42)

print(len(X_train), len(X_valid), len(X_test))


#Plot distributions of labels - make sure split preserves variation in rmsd
labels = [y_train, y_valid, y_test]
names = ['y_train', 'y_valid', 'y_test']
for i in range(0,3):
	plot_distr(labels[i], names[i], out_dir)




