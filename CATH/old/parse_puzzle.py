#! /usr/bin/env python3
# -*- coding: utf-8 -*-


import argparse
import sys
import os
import pandas as pd
import subprocess
import pdb


#Arguments for argparse module:
parser = argparse.ArgumentParser(description = '''A program that parses .dist ouput from tree-puzzle.''')
 
parser.add_argument('file_path', nargs=1, type= str,
                  default=sys.stdin, help = 'path to .dist file.')



def get_pairwise_dist(file_path):
	'''A function that gets the uids and the corresponding scores
	and returns them as lists.
	'''

	uids = [] #List with unique ids
	all_distances = [] #List with distances in order of ids in uids
	with open(file_path) as file:
		for line in file:
			line = line.rstrip()
			line = line.split(" ") #split on double space
			line = list(filter(None, line)) #Filter away empty strings

			if len(line)<2:
				number_of_seqs = int(line[0])
				continue
			else:
				for item in line:
					if '|' in item:
						uid = item.split("|")[0]
						uids.append(uid)
						new_distances = [] #Save distances for this uid
					else:
						new_distances.append(item)
				if len(new_distances) == number_of_seqs:
					all_distances.append(new_distances)
	return(uids, all_distances)



def format_distances(uids, all_distances):
	'''A function that formats all pairwise distances
	'''

	end = len(uids) #Number of uids
	count = 0 #Keep track of number of distances
	print('uid1' + '\t' + 'uid2' + '\t' + 'MLdistance')
	for i in range(0,end):
		for j in range(i+1,end):
			count+=1
			print(uids[i] + '\t' + uids[j] + '\t' + all_distances[i][j])


	number_possible_pairs = int((end/2)*(end-1))
	if number_possible_pairs != count:
		raise ValueError(str(count)+' alignmets were made out of '+ str(number_possible_pairs) + ' possible pairs.')
	return None


#Main program
args = parser.parse_args()
file_path = args.file_path[0]
(uids, all_distances)=get_pairwise_dist(file_path) #Get uids and distances
format_distances(uids, all_distances)
