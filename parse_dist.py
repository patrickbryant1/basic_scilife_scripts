#! /usr/bin/env python3
# -*- coding: utf-8 -*-


import argparse
import sys
import os
import pandas as pd
import glob
import pdb


#Arguments for argparse module:
parser = argparse.ArgumentParser(description = '''A program that parses .dist ouput from a 4x4 quartet 
												  puzzle matrix from tree-puzzle, where the ML distances from
												  duplicates of two aligned sequences have been calculated.''')
 
parser.add_argument('dir_path', nargs=1, type= str,
                  default=sys.stdin, help = '''path to directory with .dist files. The files should have names
                  on the form: uid1_uid2.phy.dist''')



def print_pairwise_dist(dir_path):
	'''A function that gets the uids and the corresponding scores
	and prints them in .tsv format.
	'''
	print('uid1' + '\t' + 'uid2' + '\t' + 'MLdistance')
	count = 0 #Keep track of number of files used
	os.chdir(dir_path)
	for file in glob.glob("*.dist"):
		#count+=1
		#print(count)

		name = file.split('.') #Split filename on .
		uids = name[0].split('_') #Separate uids
		uid_1 = uids[0] #Get uid_1
		uid_2 = uids[1] #Get uid_2

		with open(file) as file:
			for line in file:
				line = line.rstrip()
				line = line.split(" ") #split on double space
				line = list(filter(None, line)) #Filter away empty strings

				if len(line)>2:
					seq_dist = line[-1] #Get ML evolutionary distance between sequences
					print(uid_1 + '\t' + uid_2 + '\t' + seq_dist) #Print in tsv
					break

	return None



#Main program
args = parser.parse_args()
dir_path = args.dir_path[0]
print_pairwise_dist(dir_path) #Get uids and distances
