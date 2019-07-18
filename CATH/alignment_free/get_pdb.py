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

parser.add_argument('outdir', nargs=1, type= str,
                  default=sys.stdin, help = 'output directory.')
parser.add_argument('address', nargs=1, type= str,
                  default=sys.stdin, help = 'Web adress to download from.') #www.cathdb.info/version/v4_2_0/api/rest/id/

#FUNCTIONS
def get_pdb(uids, address):
	#Shuffle uids to make sure there is no selective order in comparisons within H-groups
	random.Random(2).shuffle(uids)
	end = min(len(uids), 15)
	for i in range(0,end): #Get 2-15 uid pdb files
		#Get pdb file
		subprocess.call(["wget",address+uids[i]+'.pdb'])


	return None



#####MAIN#####
args = parser.parse_args()

outdir = args.outdir[0]
address = args.address[0]

#Read df - contains H-groups with at least 2 uids filtered on 
#techinique: X-ray and resolution: <=3.5 Å
df = pd.read_csv('/home/p/pbryant/pfs/evolution/CATH/alignment_free/above2.csv')
H_group = outdir.split('/')[-2]
h_df = df[df['H_group']==H_group]
uids = [*h_df['uid']]
get_pdb(uids, address)
