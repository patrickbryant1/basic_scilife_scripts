#! /usr/bin/env python3
# -*- coding: utf-8 -*-


import argparse
import sys
import pdb

#Arguments for argparse module:
parser = argparse.ArgumentParser(description = '''A program that parses TMalign output.''')


parser.add_argument('align_file', nargs=1, type= str,
                  default=sys.stdin, help = 'path to file with TMalign output.')

#Functions
def get_pairwise_dist(file_path):
	'''A function that gets the uids and the corresponding scores
	and prints them in tsv.
	'''

	uid_pairs = [] #List with unique ids

	#df.loc[df['pdb'] == '1zmq']
	print('uid1' + '\t' + 'uid2' + '\t' + 'RMSD')
	with open(file_path) as file:
		for line in file:
			if 'Name of Chain' in line:
				line = line.rstrip() #remove \n
				line = line.split("/") #split on /
				uid = line[-1].split(".")[0] #Get uid
				
				uid_pairs.append(uid)

			if 'RMSD=' in line:
				line = line.rstrip() #remove \n
				line = line.split(",") #split on ,
				RMSD = line[1].split(' ')[-1] #Split on space

				print(uid_pairs[0] + '\t' + uid_pairs[1] + '\t' + str(RMSD))
				uid_pairs = [] #reset list of pairs



	return None


#Main program
args = parser.parse_args()

align_file = args.align_file[0]

#Get pairwise RMSD and print in tsv
get_pairwise_dist(align_file)
