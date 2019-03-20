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
parser = argparse.ArgumentParser(description = '''A program that groups uids into their h-groups.''')
 
parser.add_argument('H_group_file', nargs=1, type= str,
                  default=sys.stdin, help = 'path to file with H-groups.')

parser.add_argument('domain_file', nargs=1, type= str,
                  default=sys.stdin, help = 'path to file ECOD domain description.')


args = parser.parse_args()

H_group_file = args.H_group_file[0]
domain_file = args.domain_file[0]


def read_groups(H_group_file):
	'''Read X_group.H_group s into list
	'''

	H_groups = [] #Store H-groups
	with open(H_group_file) as file:
		for line in file:
			line = line.rstrip() #remove \n
			H_groups.append(line)

	return H_groups

#Functions
def group_ids(domain_file, H_groups):
	'''A function that gets the ids for each H-group and writes
	them to a file called X_group.H_group.txt
	'''
	count_H_groups = 0 #Count H groups
	
	for H_group in H_groups:
		uids = [] #Save uids
		pdb_ids = [] #Save pdb_ids
		count_H_groups +=1
		
		x_group = str(H_group).split('.')[0] #family level
		hom = str(H_group).split('.')[1] #homology level
		file_name = H_group +'.txt'
		with open(domain_file) as file:
			for line in file:
				line = line.rstrip() #remove \n
				
				if line[0] != '#': #Comment lines, include meta
					line = line.split("\t") #split on tab
					match_group = line[3].split('.')[0:2]
					if x_group == match_group[0]:
						if hom == match_group[1]:
							uid = line[0]
							pdb_id = line[1]
							uids.append(uid)
							pdb_ids.append(pdb_id)
				
			#After going through the entire file, the matched uids are written to a file
			write_file(file_name, uids, pdb_ids)
			print(H_group)
	
	print(count_H_groups)
	return None


def write_file(file_name, uids, pdb_ids):
	'''Write uids in same homology group to file
	'''
	try:
		with open(file_name, "w") as file:
			for i in range(0, len(uids)):
				file.write(uids[i]+ '\t' + pdb_ids[i] + '\n')

	except:
		raise IOerror('Could not write file: ' + file_name)

	return None


#####MAIN PROGRAM#####
H_groups = read_groups(H_group_file)
group_ids(domain_file, H_groups)
