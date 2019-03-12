#! /usr/bin/env python3
# -*- coding: utf-8 -*-


import argparse
import sys
import pandas as pd
import subprocess
import pdb

#Arguments for argparse module:
parser = argparse.ArgumentParser(description = '''A program that downloads pdb structures based on pdb ids.''')
 
parser.add_argument('file_path', nargs=1, type= str,
                  default=sys.stdin, help = 'path to pdb ids file. Structure: #uid	pdb_id.')

args = parser.parse_args()

file_path = args.file_path[0]

#Read tsv file as pandas dataframe
df = pd.read_csv(file_path, sep='\t')#get_ids

pdb_ids = df['pdb']

#wget http://www.rcsb.org/pdb/files/1ZMQ.pdb.gz
database_path = 'http://www.rcsb.org/pdb/files/'
for i in pdb_ids:
	subprocess.call(["wget",database_path+i+'.pdb.gz'])


