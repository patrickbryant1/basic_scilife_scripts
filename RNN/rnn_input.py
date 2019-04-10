#! /usr/bin/env python3
# -*- coding: utf-8 -*-

import numpy as np
import pandas as pd
import pdb
def read_tsv(tsv_file, threshold):
	'''Read tsv file format containing: uid1 \t uid2 \t ML distance \t RMSD distance
	'''

	uids = [] #Store uids
	rmsd_dists_t = [] #Store rmsd distances with ML_dists under threshold
	rmsd_dists = [] #Store rmsd distances
	
	with open(tsv_file) as file:
		for line in file:
			line = line.rstrip() #Remove newlines
			line = line.split("\t")
			ML_dist = round(float(line[2]), 2)
			rmsd_dists.append(float(line[3]))
			if ML_dist <= threshold:
				rmsd_dists_t.append(float(line[3]))
				
				uids.append(line[0]+'_'+line[1])				
			else:
				continue

	
	return(uids, rmsd_dists_t, rmsd_dists)


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


def make_dict(file_name, dictionary, accessibilities):
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

	
	encoding.append(len(encoding))

	return(encoding, dictionary, accessibilities)