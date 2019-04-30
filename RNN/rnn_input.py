#! /usr/bin/env python3
# -*- coding: utf-8 -*-

import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
from collections import Counter

import pdb

#ENCODING INFO

AMINO_ACIDS = { 
'A':1,
'R':2,
'N':3,
'D':4,
'C':5,
'E':6,
'Q':7,
'G':8,
'H':9,
'I':10,
'L':11,
'K':12,
'M':13,
'F':14,
'P':15,
'S':16,
'T':17,
'W':18,
'Y':19,
'V':20,
'X':21, #UNKNOWN
'-':22
}


SECONDARY_STR = {
'G':1,
'H':2,
'I':3,
'T':4,
'E':5,
'B':6,
'S':7,
'C':8,
' ':9,
'-':10
}

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
  			line = line.rstrip()
  			locations.append(line)



	return locations


def get_encodings(file_name, accessibilities, structures, letters, seqlens):
	'''Make a book of all encodings (sentences)
	'''

	aa1 = []
	str1 = []
	acc1 = []
	aa2 = []
	str2 = []
	acc2 = []

	encoding = [aa1, str1, acc1, aa2, str2, acc2] #Save int encoding of sequence

	with open(file_name) as file:
		for line in file:
			line = line.split(',')
			aa1.append(AMINO_ACIDS[line[0]])
			str1.append(SECONDARY_STR[line[1]])
			acc1.append(int(line[2]))
			aa2.append(AMINO_ACIDS[line[3]])
			str2.append(SECONDARY_STR[line[4]])
			acc2.append(int(line[5]))

			#Save all info to see distribtuions
			letters.append(line[0])
			letters.append(line[3])
			structures.append(line[1])
			structures.append(line[4])
			accessibilities.append(int(line[2]))
			accessibilities.append(int(line[5]))

			
	encoding = [aa1, str1, acc1, aa2, str2, acc2]
	#For length distributions
	seqlens.append(len(encoding[0]))

	return(encoding, accessibilities, structures, letters, seqlens)

def encoding_distributions(chart_type, encoding_info, title, bins, out_dir, name, scale):
	'''Look at distributions of encoded info
	'''
	plt.title(title)

	if chart_type == 'bar':
		counts = Counter(encoding_info)

		x = []
		y = []
		for key in counts:
			x.append(key)
			y.append(counts[key])

		
		x_pos = np.arange(len(x))
		plt.bar(x_pos , y, log = scale)

		plt.xticks(x_pos, x)


	if chart_type == 'hist':
		plt.hist(encoding_info, bins = bins, log = scale)



	plt.savefig(out_dir+name+'.png')
	plt.close()


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
