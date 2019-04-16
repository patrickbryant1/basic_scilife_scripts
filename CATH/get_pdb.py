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



def get_structures(address, uids, filter_ids, H_group):
	'''Download .pdb structure, run TMalign and check sequence identity
	'''

	downloaded_uids = [] #Keep track of ids that have been downloaded

        selected_uids = [] #uids to align
	passed_uids = [] #save uids that have passed all steps

	for i in range(0, len(uids)):
            
            if len(selected_uids) == 10:
                #Make alignment of all of these
                status = align(selected_uids)
                #If one fails, pop this and continue
                #Save the passed uids to a special list

            if uids[i][0:4].upper() in filter_ids: #Make check on pdb search
		pdb.set_trace()    
                downloaded_ids.append(uids[i])
		subprocess.call(["wget",address+selected_uids[i]+'.pdb'])
            
	    else:
                if len(selected_uids) < number_of_uids: #If all uids have not been used
		    selected_uids.insert(i+1, uids[len(selected_uids))  #remove and insert new - continue on the inserted one somehow
		else:
                    print('Not enough uids matching criteria in: ', H_group, ' Pos: ', i+1) 



        
	return None


def align(selected_uids, TMalign):
    '''Run TMalign on file pair and extract sequence identity,
    remove file2 if 90 % or above
    '''
    
    
    count = 0 #Keep track of number of alignments made
    end = len(selected_uids)

    for i in range(0, end):
            structure_i = dir_path+uids[i] #Get structure i
            for j in range(i+1, end):
                    structure_j = dir_path+uids[j] #Get structure j
                    subprocess.call([TMalign, structure_i , structure_j , '-a'])
                    count+=1
    
    return None

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

downloaded_ids = get_structures(address, uids, filter_ids, H_group)

#Print downloaded ids
print(downloaded_ids)
pdb.set_trace()



