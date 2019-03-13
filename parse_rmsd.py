#! /usr/bin/env python3
# -*- coding: utf-8 -*-


import argparse
import sys
import pandas as pd
import subprocess
import pdb

#Arguments for argparse module:
parser = argparse.ArgumentParser(description = '''A program that runs parses TMalign output.''')
 
parser.add_argument('id_file', nargs=1, type= str,
                  default=sys.stdin, help = 'path to pdb ids file. Structure: #uid	pdb_id.')

parser.add_argument('align_file', nargs=1, type= str,
                  default=sys.stdin, help = 'path to file with TMalign output.')

#Functions
def get_pairwise_dist(file_path, df):
	'''A function that gets the uids and the corresponding scores
	and returns them as lists.
	'''

	uid_pairs = [] #List with unique ids

	#df.loc[df['pdb'] == '1zmq']

	with open(file_path) as file:
		for line in file:
			if 'Name of Chain' in line:
				line = line.rstrip() #remove \n
				line = line.split("/") #split on /
				pdb_id = line[-1].split(".")[0] #Get pdb id
				uid =df.loc[df['pdb'] == pdb_id]['#uid'].values[0] #Get corresponding uid
				
				uid = str(uid).zfill(9)
				uid_pairs.append(uid)

			if 'RMSD' in line:
				line = line.rstrip() #remove \n
				line = line.split(",") #split on ,
				RMSD = line[1].split(' ')[-1] #Split on space

				print(uid_pairs[0] + '\t' + uid_pairs[1] + '\t' + str(RMSD))
				uid_pairs = [] #reset list of pairs



	return None

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
	print(str(count)+' alignmets were made out of '+ str(number_possible_pairs) + ' possible pairs.')
	return None


#Main program
args = parser.parse_args()

id_file = args.id_file[0]
align_file = args.align_file[0]

#Read tsv file as pandas dataframe
df = pd.read_csv(id_file, sep='\t')#get_ids

pdb_ids = df['pdb']
uids = df['#uid']

get_pairwise_dist(align_file, df)
