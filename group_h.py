#! /usr/bin/env python3
# -*- coding: utf-8 -*-


import argparse
import sys
import os
import pandas as pd
import subprocess
import pexpect
import pdb


#Arguments for argparse module:
parser = argparse.ArgumentParser(description = '''A program that takes aligned residue pairs from the structural
								alignment from TMalign and runs tree-puzzle on them.''')
 
parser.add_argument('H_group_file', nargs=1, type= str,
                  default=sys.stdin, help = 'path to file with H-groups.')

parser.add_argument('domain_file', nargs=1, type= str,
                  default=sys.stdin, help = 'path to file ECOD domain description.')


args = parser.parse_args()

H_group_file = args.H_group_file[0]
domain_file = args.domain_file[0]


H_group_df = pd.read_csv(H_group_file, sep='\n')#get H groups

#Functions
def group_ids(domain_file, H_group_df):
	'''A function that gets the ids for each H-group and writes
	them to a file called F_group.H_group.txt
	'''
	count_H_groups = 0 #Count H groups
	
	for H_group in H_group_df['H_group']:
		uids = [] #Save uids
		count_H_groups +=1
			
		fam = str(H_group).split('.')[0] #family level
		hom = str(H_group).split('.')[1] #homology level
		with open(domain_file) as file:
			for line in file:
				line = line.rstrip() #remove \n
				
				if line[0] != '#': #Comment lines, include meta
					line = line.split("\t") #split on tab
					match_group = line[3].split('.')[0:2]
					if fam == match_group[0]:
						if hom == match_group[1]:
							uid = line[0]
							uids.append(uid)
				
			#After going through the entire file, the matched uids are written to a file
			write_file(H_group, uids)
			print(count_H_groups)
			print(H_group)
				
	
	return None


def write_file(H_group, uids):
	'''Write uids in same homology group to file
	'''
	file_name = str(H_group) + '.txt'
	with open(file_name, "w") as file:
		for uid in uids:
			file.write(uid+'\n')

	return None


group_ids(domain_file, H_group_df)
