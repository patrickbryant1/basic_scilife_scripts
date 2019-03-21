#! /usr/bin/env python3
# -*- coding: utf-8 -*-


import argparse
import sys
import pandas as pd
import subprocess
import pdb
import glob

#Arguments for argparse module:
parser = argparse.ArgumentParser(description = '''A program that runs TMalign for pdb structures.''')
 
parser.add_argument('dir_path', nargs=1, type= str,
                  default=sys.stdin, help = 'path to base directory with .pdb structures (include / in end).')

parser.add_argument('cutoff', nargs=1, type= str,
                  default=sys.stdin, help = 'Number of structures to use for comparison')


###FUNCTIONS###

def read_ids(dir_path):
	'''Read pdb_ids into list
	'''

	pdb_ids = [] #Store pdb ids

	for text_file in glob.glob(dir_path + '*.txt'):
		with open(text_file) as file:
			for line in file:
				line = line.rstrip() #remove \n
				pdb_id = line.split('\t')[1] #Get pdb id
				pdb_ids.append(pdb_id)

	return pdb_ids

def get_pdb_file(uid):
	'''A function for formatting the path to the .pdb file
	as specified in the ECOD directory structure based on the uid.
	'''
	
	step_1 = uid[0:5]
	step_2 = uid
	step_3 = uid + '.pdbnum.pdb'

	file_path = step_1 + '/' + step_2 + '/' + step_3

	return file_path
	
def align_structures(pdb_ids, dir_path):
	'''Do structural alignment with TMalign
        '''
   
	count = 0 #Keep track of number of alignments made
	end = len(uids)

	for i in range(0, end):
		structure_i = dir_path+get_pdb_file(uids[i])
		for j in range(i+1, end):
			subprocess.call(["/home/pbryant/TMalign", structure_i , structure_j , '-a'])
			count+=1


	number_possible_pairs = int((end/2)*(end-1))
	print(str(count)+' alignmets were made out of '+ str(number_possible_pairs) + ' possible pairs.')
	
	return None

###MAIN###	
args = parser.parse_args()

dir_path = args.dir_path[0]
cutoff = args.cutoff[0]

#Get uids
pdb_ids = read_ids(dir_path)
align_structures(pdb_ids, dir_path)
