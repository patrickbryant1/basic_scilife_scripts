#! /usr/bin/env python3
# -*- coding: utf-8 -*-


import argparse
import sys
import subprocess
import pdb
import glob

#Arguments for argparse module:
parser = argparse.ArgumentParser(description = '''A program that runs TMalign for pdb structures.''')
 
parser.add_argument('dir_path', nargs=1, type= str,
                  default=sys.stdin, help = 'path to base directory with .pdb structures (include / in end).')

parser.add_argument('TMalign_path', nargs=1, type= str,
                  default=sys.stdin, help = 'path to TMalign.')


###FUNCTIONS###

def read_ids(dir_path):
	'''Read pdb_ids into list
	'''

	uids = [] #Store uids

	for file_name in glob.glob(dir_path + '*.pdb'):
		
		uid = file_name.split('/')[-1] #The last part of the path is the filename
		uids.append(uid) #append to uids

	return uids

	
def align_structures(uids, dir_path, TMalign):
	'''Do structural alignment with TMalign
    '''
   
	count = 0 #Keep track of number of alignments made
	end = len(uids)

	for i in range(0, end):
		structure_i = dir_path+uids[i] #Get structure i
		for j in range(i+1, end):
			structure_j = dir_path+uids[j] #Get structure j
			subprocess.call([TMalign, structure_i , structure_j , '-a'])
			count+=1


	number_possible_pairs = int((end/2)*(end-1))
	print(str(count)+' alignmets were made out of '+ str(number_possible_pairs) + ' possible pairs.')
	
	return None

###MAIN###	
args = parser.parse_args()

dir_path = args.dir_path[0]
TMalign = args.TMalign_path[0]


#Get uids
uids = read_ids(dir_path)
#Align structures
align_structures(uids, dir_path, TMalign)
