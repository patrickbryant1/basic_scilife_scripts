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
parser.add_argument('filter_file', nargs=1, type= str,
                  default=sys.stdin, help = 'path to file that contains newline separated pdb ids from a pdb search.')
parser.add_argument('address', nargs=1, type= str,
                  default=sys.stdin, help = 'Web adress to download from.') #www.cathdb.info/version/v4_2_0/api/rest/id/
parser.add_argument('output_dir', nargs=1, type= str,
                  default=sys.stdin, help = 'output directory.')


#FUNCTIONS
def read_newline(file_name):
	'''Read newline separated file contents into list
	'''

	
	contents = [] #Store contents

	with open(file_name) as file:
			
		for line in file:
			line = line.rstrip() #remove \n
			contents.append(line) #Add uid
	

	return(contents)



def get_structures(address, uids, filter_ids):
	'''Download .pdb structure
	'''

	downloaded_ids = [] #Keep track of ids that have been downloaded
	
	number_of_uids = len(uids)

	selected_uids = uids[0:10]

	for i in range(0, len(selected_uids)):

		if selected_uids[i].upper() in filter_ids:
			found = glob.glob(output_dir + uid + '*')

			if found
			downloaded_ids.append(uid)
			subprocess.call(["wget",address+uid+'.pdb'])

		else:
			selected_uids. #remove and insert new - continue on the inserted one somehow
			

	return downloaded_ids

#####MAIN#####
args = parser.parse_args()

uid_file = args.uid_file[0]
filter_file = args.filter_file[0]
address = args.address[0]
output_dir = args.output_dir[0]


#Get selected uids:
H_group = uid_file.split('/')[-1] #Get H-group (last part of path)
uids = read_newline(uid_file)

#Get pdb ids to filter on
filter_ids = read_newline(filter_file)

pdb.set_trace()

downloaded_ids = get_structures(address, uids, H_groups)

#Print downloaded ids
print(downloaded_ids)
pdb.set_trace()



