#! /usr/bin/env python3
# -*- coding: utf-8 -*-


import argparse
import sys
import pandas as pd
import subprocess
import pdb

#Arguments for argparse module:
parser = argparse.ArgumentParser(description = '''A program that runs TMalign for pdb structures.''')
 
parser.add_argument('file_path', nargs=1, type= str,
                  default=sys.stdin, help = 'path to pdb ids file. Structure: #uid	pdb_id.')

parser.add_argument('dir_path', nargs=1, type= str,
                  default=sys.stdin, help = 'path to directory with .pdb structures (include / in end.')

args = parser.parse_args()

file_path = args.file_path[0]
dir_path = args.dir_path[0]

#Read tsv file as pandas dataframe
df = pd.read_csv(file_path, sep='\t')#get_ids

pdb_ids = df['pdb']

#Do structural alignment with TMalign
count = 0
position = 0
while len(pdb_ids)>1: #While structures are left to be aligned
	for i in range(position+1, len(pdb_ids)):
		subprocess.call(["/home/pbryant/TMalign", dir_path+pdb_ids[position]+'.pdb', dir_path+pdb_ids[i]+'.pdb' ])
		count+=1
		print(count)
	print(pdb_ids)
	pdb_ids.pop(position)
	position+=1
	

