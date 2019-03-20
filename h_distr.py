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
									uids in each H-group in ECOD.''')
 
parser.add_argument('dir_path', nargs=1, type= str,
                  default=sys.stdin, help = 'path to file with H-groups.')


args = parser.parse_args()

dir_path= args.dir_path[0]



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


#######################MAIN################################
(h_counts, h_group, id_count) = compute_stats(dir_path)

log_id_count = numpy.log10(id_count)

plot_hist(log_id_count, 100)
pdb.set_trace()


