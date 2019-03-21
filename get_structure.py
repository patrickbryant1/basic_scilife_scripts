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
                  default=sys.stdin, help = 'path to pdb ids file directory. Structure of files: #uid	pdb_id.')
parser.add_argument('address', nargs=1, type= str,
                  default=sys.stdin, help = 'Web adress to download from.')
parser.add_argument('output_dir', nargs=1, type= str,
                  default=sys.stdin, help = 'output directory.')
parser.add_argument('selected_files', nargs=1, type= str,
                  default=sys.stdin, help = 'Newline separated file with file names for files with ids to get structures for.')
parser.add_argument('n_entries', nargs=1, type= int,
                  default=sys.stdin, help = 'Number of entreis to get for each H-group.')




def read_selected(selected_files):
	'''Read newline separated filenames into list
	'''

	selected = []
	with open(selected_files) as file:
		for line in file:
			line = line.rstrip() #remove \n
			selected.append(line)

	return(selected)

#Read tsv file 
def read_groups(file_path):
	'''Read tsv file into list
	'''

	uids = [] #Store uids
	pdb_ids = [] #Store pdb_ids

	with open(file_path) as file:
		for line in file:
			line = line.rstrip() #remove \n
			line = line.split('\t')
			uids.append(line[0])
			pdb_ids.append(line[1])

	return(uids, pdb_ids)

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
selected_files = args.selected_files[0]
n_entries = args.n_entries[0]

#Get selected file names:
selected = read_selected(selected_files)

downloaded_ids = [] #Keep track of ids that have been downloaded
for file_name in selected:
	downloaded_ids = get_structures(address, file_name, input_dir, output_dir, n_entries, downloaded_ids)

#Print downloaded ids
print(downloaded_ids)
for downloaded_id in downloaded_ids:
	print(downloaded_id)



