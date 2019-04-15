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
parser = argparse.ArgumentParser(description = '''A program that downloads pdb structures based on CATH uids (domain ids)
												from H-groups that have at least 10 entries.
												It then runs TMalign on all pairs of 10 randomly selected entries. If a paired
												alignments should happen to have above 90% sequence identity, the second uid
												is dropped and a new domain structure downloaded, if available. Otherwise
												the whole H-group is dropped''')
 
parser.add_argument('input_dir', nargs=1, type= str,
                  default=sys.stdin, help = 'path to CATH ids file directory.')
parser.add_argument('address', nargs=1, type= str,
                  default=sys.stdin, help = 'Web adress to download from.')
parser.add_argument('output_dir', nargs=1, type= str,
                  default=sys.stdin, help = 'output directory.')

#www.cathdb.info/version/v4_2_0/api/rest/id/


#FUNCTIONS

#Functions
def read_tsv(file_path):
	'''Read ids and H-groups into lists 
	'''

	uids = [] #Store ids
	H_groups = [] #Store H-groups

	with open(file_path) as file:
		for line in file:
			line = line.rstrip() #remove \n
			line = line.split(',')
			uid = line[0]
			H_group = line[1:]

			#Add id into right h_group
			if H_groups:
				i = 0 #Reset index
				found = False #Keep track of matches
				while i < len(H_groups):
					#If a match is found
					if H_groups[i] == H_group:
						uids[i].append(uid)
						found = True
						print(len(H_groups))
						break #Break out of loop
					i+=1

				#If no match is found
				if found == False:
					H_groups.append(H_group)
					uids.append([uid])
			else:
				H_groups.append(H_group)
				uids.append([uid])


	return uids, H_groups

def select_n_random(uids, n):
	'''Select n random uids from each H_group
	'''

	selected = random.sample(uids, n)

	return selected

for i in range(0, len(uids)):
	uid_counts.append(len(uids[i]))
	if len(uids[i])>=n:
		over_n.append(H_groups[i])
		selected = select_n_random(uids[i], n)
		n_random.append(selected)

def get_structures(address, uids, H_groups):
	'''Download .pdb structure
	'''

	downloaded_ids = [] #Keep track of ids that have been downloaded
	
	for i in range(0, len(H_groups)):
		dir_name = H_groups[i] #Get directory name (C.A.T.H)
		subprocess.call(['mkdir', dir_name])

		for uid in uids[i]:
			downloaded_ids.append(uid)
			subprocess.call(["wget",address+uid+'.pdb'])

		for file in glob.glob(output_dir + '*.pdb'):
			shutil.move(file, output_dir+dir_name)
			

	return downloaded_ids

#####MAIN#####
args = parser.parse_args()

input_dir = args.input_dir[0]
address = args.address[0]
output_dir = args.output_dir[0]


#Get selected file names:
(uids, H_groups) = read_ids(input_dir)

#Make check 
if len(uids) != len(H_groups):
	raise ValueError('There are not uids for every H-group!')


downloaded_ids = get_structures(address, uids, H_groups)

#Print downloaded ids
print(downloaded_ids)
pdb.set_trace()



