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
from collections import Counter


#Arguments for argparse module:
parser = argparse.ArgumentParser(description = '''A program that investigates the distribution of
									 each H-group in CATH.''')
 
parser.add_argument('file_path', nargs=1, type= str,
                  default=sys.stdin, help = 'path to file with H-groups.')


args = parser.parse_args()

file_path= args.file_path[0]


#Functions
def read_tsv(file_path):
	'''Read ids and H-groups into list
	'''

	uids = [] #Store ids
	H_groups = [] #Store H-groups

	with open(file_path) as file:
		for line in file:
			line = line.rstrip() #remove \n
			line = line.split(',')
			uid = line[0]
			H_group = line[1]

			uids.append(uid)
			H_groups.append(H_group)

	return uids, H_groups


def plot_hist(id_count, bins):
	plt.hist(id_count, normed = True, bins = bins)
	plt.show()

	return None

def select_n(count_h_list, unique_h_list, n):
	'''Select H-groups with at least n entries.
	'''

	over_n = [] #Store H_groups with over n entries
	over_n_count = [] #Store counts for H_groups with over n entries
	for i in range(0, len(count_h_list)):
		if count_h_list[i] >= n:
			over_n.append(unique_h_list[i])
			over_n_count.append(count_h_list[i])
			print(unique_h_list[i])

	print('Number of H-groups with at least ' + str(n) + 'entries: ' + str(len(over_n)))

	return(over_n, over_n_count)


#######################MAIN################################
uids, H_groups = read_tsv(file_path)

count_h = Counter(H_groups).values() #Count occurrence of each h group
unique_h = Counter(H_groups).keys() #Get all unique h groups

#Turn dicts into lists
count_h_list = []
unique_h_list = []
for i in count_h:
	count_h_list.append(i)

for j in unique_h:
	unique_h_list.append(j)


log_h_count = numpy.log10(count_h_list)
(over_n, over_n_count) = select_n(count_h_list, unique_h_list, 10)

pdb.set_trace()


