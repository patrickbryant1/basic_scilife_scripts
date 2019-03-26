#! /usr/bin/env python3
# -*- coding: utf-8 -*-


import argparse
import sys
import subprocess
import pdb
import glob
import shutil
import random


#Arguments for argparse module:
parser = argparse.ArgumentParser(description = '''A program that downloads pdb structures based on pdb ids.''')
 
parser.add_argument('input_dir', nargs=1, type= str,
                  default=sys.stdin, help = 'path to CATH ids file directory.')
parser.add_argument('address', nargs=1, type= str,
                  default=sys.stdin, help = 'Web adress to download from.')
parser.add_argument('output_dir', nargs=1, type= str,
                  default=sys.stdin, help = 'output directory.')



#FUNCTIONS

def read_ids(input_dir):
	'''Read newline separated fileS WITH CATH ids into list
	'''

	pdb_ids = []

	for file_name in glob.glob(input_dir + '*'):
		with open(file_name) as file:
			for line in file:
				line = line.rstrip() #remove \n
				pdb_id = line[0:4] #The first four chars are the pdb id 
				if pdb_id not in pdb_ids:
					print(pdb_id)
					pdb_ids.append(pdb_id)

	return(pdb_ids)

def write_ids(destination, file_name, positions, uids, pdb_ids):
	'''Write uids and pdb_ids to .txt file
	'''

	with open(destination+'/'+file_name, "w") as file:
		for i in positions:
			file.write(uids[i]+ '\t' + pdb_ids[i] + '\n')
	

	return None

def get_structures(address, file_name, input_dir, output_dir, n_entries, downloaded_ids):
	'''Download .pdb structure and group into directory
	'''
	file_path = input_dir + file_name
	(uids, pdb_ids) = read_groups(file_path)
	positions = random.sample(range(0, len(uids)), n_entries) #Get n_entries random positions to download files for these ids

	dir_name = file_name.split('.txt')[0] #Get directory name (X.H group)
	subprocess.call(['mkdir', dir_name])

	for i in positions:
		downloaded_ids.append(pdb_ids[i])
		subprocess.call(["wget",address+pdb_ids[i]])

	for file in glob.glob(output_dir + 'str*'):
		shutil.move(file, output_dir+dir_name)

	write_ids(output_dir+dir_name, file_name, positions, uids, pdb_ids)	

	
	return downloaded_ids

#####MAIN#####
args = parser.parse_args()

input_dir = args.input_dir[0]
address = args.address[0]
output_dir = args.output_dir[0]


#Get selected file names:
pdb_ids = read_ids(input_dir)
pdb.set_trace()

downloaded_ids = [] #Keep track of ids that have been downloaded
for file_name in selected:
	downloaded_ids = get_structures(address, file_name, input_dir, output_dir, n_entries, downloaded_ids)

#Print downloaded ids
print(downloaded_ids)
for downloaded_id in downloaded_ids:
	print(downloaded_id)



