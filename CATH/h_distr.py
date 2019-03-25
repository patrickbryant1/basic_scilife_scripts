#! /usr/bin/env python3
# -*- coding: utf-8 -*-


import argparse
import sys
import os
import pandas as pd
import subprocess
import glob
import matplotlib.pyplot as plt
import pdb
import numpy


#Arguments for argparse module:
parser = argparse.ArgumentParser(description = '''A program that investigates the distribution of
									 each H-group in CATH.''')
 
parser.add_argument('dir_path', nargs=1, type= str,
                  default=sys.stdin, help = 'path to file with H-groups.')


args = parser.parse_args()

dir_path= args.dir_path[0]


def read_ids(dir_path):
	'''Read ids and H-groups into list
	'''

	ids = [] #Store ids
	H_groups = [] #Store H-groups

	for text_file in glob.glob(dir_path + '*.txt'):
		with open(text_file) as file:
			for line in file:
				line = line.rstrip() #remove \n
				pdb_id = line.split('\t')[1] #Get pdb id
				pdb_ids.append(pdb_id)

	return pdb_ids
#Functions
def compute_stats(dir_path):
	'''A function that counts the ids for each H-group and returns
	statistics regarding them.	
	'''

	h_counts = {} #dir to save id counts of each H-group
	h_group = [] #Store h_group
	id_count = [] #Store id counts for histogram

	os.chdir(dir_path)
	for file in glob.glob("*.txt"):
		if os.stat(file).st_size == 0:
			print('empty_file: '+file)
		else:
			id_df = pd.read_csv(file, header = None, sep='\n')#get ids for H group, don't read header (no header)
			num_ids = len(id_df)
			h_counts[file] = num_ids
			h_group.append(file)
			id_count.append(num_ids)



	return(h_counts, h_group, id_count)

def plot_hist(id_count, bins):
	plt.hist(id_count, normed = True, bins = bins)
	plt.show()

	return None

def select_n(h_counts, n):
	'''Select H-groups with at least n entries.
	'''

	over_n = {}
	for key in h_counts:
		n_entries = h_counts[key]
		if n_entries >= n:
			print(key)
			over_n[key] = n_entries

	print('Number of H-groups with at least ' + str(n) + 'entries: ' + str(len(over_n)))

	return over_n


#######################MAIN################################
(h_counts, h_group, id_count) = compute_stats(dir_path)

log_id_count = numpy.log10(id_count)
over_n = select_n(h_counts, 10)

pdb.set_trace()


