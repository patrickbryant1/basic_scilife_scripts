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
									 each H-group in CATH. And writes all uids for each H-group with at least x entries 
									 into files(newline separated) in the output directory''')
 
parser.add_argument('file_path', nargs=1, type= str,
                  default=sys.stdin, help = 'Path to file with H-groups.')

parser.add_argument('outdir_path', nargs=1, type= str,
                  default=sys.stdin, help = 'Path to output directory.')

parser.add_argument('x', nargs=1, type= int,
                  default=sys.stdin, help = 'Minimum number of entries per H-group (integer).')


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
			H_group = line[1:]

			#Add id into right h_group
			if H_groups:
				i = 0 #Reset index
				found = False #Keep track of matches
				while i < len(H_groups):
					#If a match is found
					if H_groups[i] == H_group:
						uids[i].append(uid)
						found = True
						break #Break out of loop
					i+=1

				#If no match is found
				if found == False:
					H_groups.append(H_group)
					uids.append([uid])
					print(len(H_groups))
			else:
				H_groups.append(H_group)
				uids.append([uid])



	return uids, H_groups


def plot_hist(id_count, bins, title):
	plt.title(title)
	plt.xlabel('Number of uids in H-group')
	plt.ylabel('log(count) of number of H-groups with x uids')
	plt.hist(id_count, bins = bins)
	plt.yscale("log")
	plt.xscale("log")
	plt.show()


	return None



def write_selected(over_n_H, over_n_uids, outdir_path):
	'''Write the selected uids into a file named H_group
	'''

	for i in range(0, len(over_n_H)):
		#Open a file named "H-group to write to"
		H_group = over_n_H[i]
		name = str(H_group[0])+'.'+str(H_group[1])+'.'+str(H_group[2])+'.'+str(H_group[3])
		with open(outdir_path+name, 'w') as file:
			#Write uids to file
			for uid in over_n_uids[i]:
				file.write(uid+'\n')

	return None

#######################MAIN################################
args = parser.parse_args()

file_path = args.file_path[0]
outdir_path = args.outdir_path[0]
x = args.x[0]

uids, H_groups = read_tsv(file_path)

#Count uids in each H_group:
uid_counts_all = [] #Store uid counts
uid_counts_over = [] #Store uid counts over n
over_n_H = [] #Store H_groups with uid counts over n
over_n_uids = [] #Store uids with over n entries
n = x #Cutoff

for i in range(0, len(uids)):
	uid_counts_all.append(len(uids[i]))
	if len(uids[i])>=n:
		over_n_H.append(H_groups[i])
		selected = uids[i]
		over_n_uids.append(selected)
		uid_counts_over.append(len(uids[i]))

pdb.set_trace()

write_selected(over_n_H, over_n_uids, outdir_path)






