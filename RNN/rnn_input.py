#! /usr/bin/env python3
# -*- coding: utf-8 -*-

import numpy as np
import pandas as pd
import matplotlib.pyplot as plt

import pdb



#Functions for reading, exploring and visualizing input data

def read_labels(tsv_file):
	'''Read tsv file format containing: uid1 \t uid2 \t ML distance \t RMSD distance
	'''

	distance_dict = {} #Save information for each uid pair
	
	with open(tsv_file) as file:
		for line in file:
			line = line.rstrip() #Remove newlines
			line = line.split("\t") #Split on tab
			uid_pair = (line[0]+'_'+line[1]) #Get uid pair
			
			ML_dist = float(line[2]) #aa evolutionary distance from tree-puzzle
			rmsd_dist = float(line[3]) #rmsd from TMalign

			distance_dict[uid_pair] = [ML_dist, rmsd_dist]
			
				

	return(distance_dict)


def rmsd_hot(rmsd_dists):
	'''Normalize rmsd distances and convert to one-hot encodings
	'''

	#convert to numpy array
	rmsd_dists = np.array(rmsd_dists)
	#Normalize due to max rmsd
	rmsd_dists = rmsd_dists/max(rmsd_dists)
	#Multiply with 100
	rmsd_dists = rmsd_dists*100
	#Round to 0 decimal pints
	rmsd_dists = np.around(rmsd_dists, decimals = 0)
	#Convert to ints
	rmsd_dists = rmsd_dists.astype(int)
	#Convert to one-hot encoding
	rmsd_dists_hot = np.eye(101)[rmsd_dists]


	return rmsd_dists_hot

def get_locations(encode_locations):
	'''Get all file locations for encodings
	'''

	locations = []
	with open(encode_locations) as file:
  		for line in file:
  			if '/home/pbryant' in line:
  				directory = line.rstrip()[:-1]
  			else:
  				locations.append(directory+'/'+line.rstrip())



	return locations


def make_dict(file_name, dictionary, accessibilities, seqlens):
	'''Make a dictionary of all encodings
	'''

	encoding = [] #Save int encoding of sequence
	
	with open(file_name) as file:
		for line in file:
			line = line.split(',')
			acc1 = int(line[2])
			acc2 = int(line[5])
			accessibilities.append(acc1)
			accessibilities.append(acc2)

			#Sort acc to categorial (1 if above 10)
			if acc1 > 10:
				line[2] = '1'
			else:
				line[2] = '0'

			if acc2 > 10:
				line[5] = '1'
			else:
				line[5] = '0'


			word = ''.join(line[0:6])
			

			if word in dictionary:
				encoding.append(dictionary[word])
			else:
				dictionary[word] = int(len(dictionary))
				encoding.append(dictionary[word])

	seqlens.append(len(encoding))
	#encoding.append(len(encoding))

	return(encoding, dictionary, accessibilities, seqlens)


def get_labels(encodings, distance_dict):
	'''Get corresponding labels for encodings
	'''

	uids = [] #save uids
	encoding_list = [] #save encodings
	rmsd_dists = []
	ML_dists = []


	for key in encodings:
		uids.append(key)

		encoding_list.append(encodings[key])
		[ML_dist, rmsd_dist] = distance_dict[key]

		rmsd_dists.append(rmsd_dist)
		ML_dists.append(ML_dist)



	return (uids, encoding_list, rmsd_dists, ML_dists)


def word_distributions(encoding_list, bins, out_dir, name):
	'''Count the occurence of all words in all encodings in encoding_list
	'''

	cat_encodings = [j for i in encoding_list for j in i]
		
	plt.hist(cat_encodings, bins = bins, log = True)
	plt.savefig(out_dir+name+'_word_distr.png')
	plt.close()

	return None


def label_distr(y, name, out_dir):
	'''plot distribution of labels (normalized rmsd values)
	'''

	#Convert back to ints
	y_ints = np.argmax(y, axis = 1)

	plt.hist(y_ints, bins = 101)
	plt.savefig(out_dir+name+'_label_distr.png')
	plt.close()
	
	#plt.show()
	return None
