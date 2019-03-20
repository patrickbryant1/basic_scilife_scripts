#! /usr/bin/env python3
# -*- coding: utf-8 -*-


import argparse
import sys
import subprocess
import pdb

#Arguments for argparse module:
parser = argparse.ArgumentParser(description = '''A program that downloads pdb structures based on pdb ids.''')
 
parser.add_argument('file_path', nargs=1, type= str,
                  default=sys.stdin, help = 'path to pdb ids file. Structure: #uid	pdb_id.')
parser.add_argument('address', nargs=1, type= str,
                  default=sys.stdin, help = 'Web adress to download from.')


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

#####MAIN#####
args = parser.parse_args()

file_path = args.file_path[0]
address = args.address[0]

(uids, pdb_ids) = read_groups(file_path)

#wget http://www.rcsb.org/pdb/files/1ZMQ.pdb.gz
database_path = 'http://www.rcsb.org/pdb/files/'
for i in pdb_ids:
	subprocess.call(["wget",database_path+i+'.pdb.gz'])


