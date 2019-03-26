#! /usr/bin/env python3
# -*- coding: utf-8 -*-


import argparse
import sys
import os
import matplotlib.pyplot as plt
import numpy
import random

import pdb

#Arguments for argparse module:
parser = argparse.ArgumentParser(description = '''A program that investigates the distribution of
									 each H-group in CATH. And writes n randomly selected entries 
									 for each H-group into files(newline separated) in the output directory''')
 
parser.add_argument('file_path', nargs=1, type= str,
                  default=sys.stdin, help = 'Path to file with H-groups.')

parser.add_argument('outdir_path', nargs=1, type= str,
                  default=sys.stdin, help = 'Path to output directory.')


#Functions
def read_tsv(file_path):
	'''Read ids and H-groups into lists 
	'''

	uids = [] #Store ids
	H_groups = [] #Store H-groups

	with open(file_path) as file:
		for line in file:
			line = line.rstrip() #remove \n
			line = line.split(',')
			uid = line[0]
			H_group = line[1]

			#Add id into right h_group
			if H_groups:
				i = 0 #Reset index
				found = False #Keep track of matches
				while i < len(H_groups):
					#If a match is found
					if H_groups[i] == H_group:
						uids[i].append(uid)
						found = True
						break
					i+=1

				#If no match is found
				if found == False:
					H_groups.append(H_group)
					uids.append([uid])
			else:
				H_groups.append(H_group)
				uids.append([uid])


	return uids, H_groups


def plot_hist(id_count, bins):
	plt.hist(id_count, normed = True, bins = bins)
	plt.show()

	return None

def select_n_random(uids, n):
	'''Select n random uids from each H_group
	'''

	selected = random.sample(uids, n)

	return selected

def write_selected(over_n, n_random, outdir_path):
	'''Write the selected uids into a file named H_group
	'''

	for i in range(0, len(over_n)):
		#Open a file named "H-group to write to"
		with open(outdir_path+str(over_n[i]), 'w') as file:
			#Write uids to file
			for uid in n_random[i]:
				file.write(uid+'\n')

	return None

#######################MAIN################################
args = parser.parse_args()

file_path = args.file_path[0]
outdir_path = args.outdir_path[0]

uids, H_groups = read_tsv(file_path)

#Count uids in each H_group:
uid_counts = [] #Store uid_counts
over_n = [] #Store H_groups with uid counts over n
n = 10 #Cutoff
n_random = [] #Store n randomly selected uids from the H_groups with over n entries
for i in range(0, len(uids)):
	uid_counts.append(len(uids[i]))
	if len(uids[i])>=n:
		over_n.append(H_groups[i])
		selected = select_n_random(uids[i], n)
		n_random.append(selected)


write_selected(over_n, n_random, outdir_path)

#Convert to log
#log_h_count = numpy.log10(count_h_list)




