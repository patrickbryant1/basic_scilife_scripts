#! /usr/bin/env python3
# -*- coding: utf-8 -*-


import argparse
import sys
import subprocess
import pdb
import shutil
import random


#Arguments for argparse module:
parser = argparse.ArgumentParser(description = '''A program that downloads pdb structures based on CATH uids (domain ids)
												from H-groups that have at least 10 entries.
												It then runs TMalign on all pairs of 10 randomly selected entries. If a paired
												alignments should happen to have above 90% sequence identity, the second uid
												is dropped and a new domain structure downloaded, if available. Otherwise
												the whole H-group is dropped''')
 
parser.add_argument('uid_file', nargs=1, type= str,
                  default=sys.stdin, help = 'path to CATH ids file.')
parser.add_argument('address', nargs=1, type= str,
                  default=sys.stdin, help = 'Web adress to download from.')
parser.add_argument('output_dir', nargs=1, type= str,
                  default=sys.stdin, help = 'output directory.')

#www.cathdb.info/version/v4_2_0/api/rest/id/


#FUNCTIONS

def read_ids(input_dir):
	'''Read newline separated fileS WITH CATH ids into list
	'''

	
	uids = [] #Store uids

	H_group = file_name.split('/')[-1] #Get H-group (last part of path)

		with open(file_name) as file:
			
			for line in file:
				uid = line.rstrip() #remove \n
				uids.append(uid) #Add uid
	

	return(uids, H_group)

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

uid_file = args.uid_file[0]
address = args.address[0]
output_dir = args.output_dir[0]


#Get selected file names:
(uids, H_group) = read_ids(uid_file)

#Make check 
if len(uids) != len(H_groups):
	raise ValueError('There are not uids for every H-group!')


downloaded_ids = get_structures(address, uids, H_groups)

#Print downloaded ids
print(downloaded_ids)
pdb.set_trace()



