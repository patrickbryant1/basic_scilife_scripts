#! /usr/bin/env python3
# -*- coding: utf-8 -*-


import argparse
import sys
import subprocess
import random
import os
import glob
import pandas as pd
import pdb



#Arguments for argparse module:
parser = argparse.ArgumentParser(description = '''Gets pdb files from CATH according to uid''')

parser.add_argument('filter_file', nargs=1, type= str,
                  default=sys.stdin, help = 'path to file that contains newline separated pdb ids from a pdb search.')
parser.add_argument('outdir', nargs=1, type= str,
                  default=sys.stdin, help = 'output directory.')
parser.add_argument('address', nargs=1, type= str,
                  default=sys.stdin, help = 'Web adress to download from.') #www.cathdb.info/version/v4_2_0/api/rest/id/

#FUNCTIONS
def read_newline(file_name):
	'''Read newline separated file contents into list
	'''
	contents = [] #Store contents

	with open(file_name) as file:

		for line in file:
			line = line.rstrip() #remove \n
			contents.append(line) #Add

	return(contents)

def get_pdb(uids, address):
	pdb.set_trace()
	#Shuffle uids to make sure there is no selective order in comparisons within H-groups
	random.Random(2).shuffle(uids)
	pdb.set_trace()
	for uid in uids:
		#Get pdb file
		subprocess.call(["wget",address+uid+'.pdb'])


	return None



#####MAIN#####
args = parser.parse_args()

outdir = args.outdir[0]
filter_file = args.filter_file[0]
address = args.address[0]

#Get filter ids
filter_ids = read_newline(filter_file)
#Read df
df = pd.read_csv('/home/p/pbryant/pfs/evolution/CATH/alignment_free/above2.csv')
H_group = outdir.split('/')[-2]
h_df = df[df['H_group']==H_group]
uids = [*h_df['uid']]
get_pdb(uids, address)
