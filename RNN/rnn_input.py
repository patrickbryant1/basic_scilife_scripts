#! /usr/bin/env python3
# -*- coding: utf-8 -*-

import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import collections
from scipy import stats

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
' ':8,
'-':9
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
			Chain1 = int(line[4]) 
			Chain2 = int(line[5])
			Aligned = int(line[6])/min(Chain1, Chain2)
			Identity = float(line[7])
			

			distance_dict[uid_pair] = [ML_dist, rmsd_dist, Chain1, Chain2, Aligned, Identity]
			
				

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


def get_encodings(file_name):
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

			

			
	encoding = [aa1, str1, acc1, aa2, str2, acc2]

	return(encoding)

def encoding_distributions(chart_type, encoding_info, title, xlabel, ylabel, bins, out_dir, name, scale, xticks):
	'''Look at distributions of encoded info
	'''
	plt.title(title)
	plt.xlabel(xlabel)
	plt.ylabel(ylabel)

	if chart_type == 'bar':
		counts = collections.Counter(encoding_info)
		#sort dicttionary
		od = collections.OrderedDict(sorted(counts.items()))
		x = []
		y = []
		for key in od:
			x.append(key)
			y.append(counts[key])

		plt.bar(x , y, log = scale)

		plt.xticks(x, xticks)


	if chart_type == 'hist':
		plt.hist(encoding_info, bins = bins, log = scale)
		mean = sum(encoding_info)/len(encoding_info)
		plt.axvline(mean, color='k', linestyle='dashed', linewidth=1)

	plt.savefig(out_dir+name+'.png')
	plt.close()


def get_labels(encodings, distance_dict, threshold):
	'''Get corresponding labels for encodings
	'''

	uids = [] #save uids
	encoding_list = [] #save encodings
	rmsd_dists = []
	ML_dists = []
	Chains = []
	Align_lens = []
	Identities = []

	letters = []
	structures = []
	accessibilities = []
	seqlens = []
	
	for key in encodings:
		uids.append(key)

		
		[ML_dist, rmsd_dist, Chain1, Chain2, Aligned, Identity] = distance_dict[key]

		if ML_dist <= threshold:
			#Save all info to see distribtuions and to feed later
			#From tree-puzzle and TMalign
			rmsd_dists.append(rmsd_dist)
			ML_dists.append(ML_dist)
			Chains.append(Chain1)
			Chains.append(Chain2)
			Align_lens.append(Aligned)
			Identities.append(Identity)

			#From encoding
			encoding = encodings[key]
			for pos in range(0, len(encoding[0])):

				letters.append(encoding[0][pos])
				letters.append(encoding[3][pos])
				structures.append(encoding[1][pos])
				structures.append(encoding[4][pos])
				accessibilities.append(encoding[2][pos])
				accessibilities.append(encoding[5][pos])

			seqlens.append(len(encoding[0]))
			encoding_list.append(encoding)

	return (uids, encoding_list, rmsd_dists, ML_dists, Chains, Align_lens, Identities, letters, structures, accessibilities, seqlens)


def word_distributions(encoding_list, bins, out_dir, name):
	'''Count the occurence of all words in all encodings in encoding_list
	'''

	cat_encodings = [j for i in encoding_list for j in i]
		
	plt.hist(cat_encodings, bins = bins, log = True)
	plt.savefig(out_dir+name+'_word_distr.png')
	plt.close()

	return None


def label_distr(x, y, name, out_dir, xlabel, ylabel):
	'''plot distribution of labels (normalized rmsd values)
	'''
	xedges, yedges = np.linspace(0, 9, 10*9), np.linspace(0, 8, 10*8)
	hist, xedges, yedges = np.histogram2d(x, y, (xedges, yedges))

	xidx = np.clip(np.digitize(x, xedges), 0, hist.shape[0]-1)
	yidx = np.clip(np.digitize(y, yedges), 0, hist.shape[1]-1)
	c = hist[xidx, yidx]
	plt.scatter(x, y, c=c)

	#Calculate line of best fit
	(slope, intercept, r_value, p_value, std_err) = stats.linregress(x, y)

	#Desciption
	plt.title('Sequence vs Structural distance' + '\n' + 'R-squared: ' + str((r_value**2).round(3)) +'|' + 'Slope: ' + str(slope.round(3)))
	plt.xlabel(xlabel)
	plt.ylabel(ylabel)
	#Line of best fit
	plt.plot(x, intercept + slope*np.array(x), 'r')
	#Colorbar
	plt.colorbar()

	plt.savefig(out_dir+name+'.png')
	plt.close()
	
	
	return None
