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
parser.add_argument('TMalign_path', nargs=1, type= str,
                  default=sys.stdin, help = 'path to TMalign.')

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



def get_structures(address, uids, filter_ids, H_group, TMalign):
	'''Download .pdb structure, run TMalign and check sequence identity
	'''

	downloaded_uids = [] #Keep track of ids that have been downloaded	
	passed_uids = [] #save uids that have passed all steps
	selected_uids = []
	for i in range(0, len(uids)):
		if len(selected_uids) == 10:
			#Make alignment of all of these
			pdb.set_trace()
			status = align(selected_uids, TMalign)
			#If one fails, pop this and continue
			#Save the passed uids to a special list

		if uids[i][0:4].upper() in filter_ids: #Make check on pdb search
		
			downloaded_uids.append(uids[i])
			selected_uids.append(uids[i])
			print(uids[i])
			subprocess.call(["wget",address+uids[i]+'.pdb'])
            
		


        
	return None


def align(selected_uids, TMalign):
    '''Run TMalign on file pair and extract sequence identity,
    remove file2 if 90 % or above
    '''


    
    
    count = 0 #Keep track of number of alignments made
    end = len(selected_uids)

    for i in range(0, end):
            structure_i ='/home/p/pbryant/pfs/evolution/CATH/'+selected_uids[i]+'.pdb' #Get structure i
            for j in range(i+1, end):
                    structure_j ='/home/p/pbryant/pfs/evolution/CATH/'+selected_uids[j]+'.pdb' #Get structure j
                    subprocess.call([TMalign, structure_i , structure_j , '-a'])
                    count+=1
    


	

    #Parse
    pdb.set_trace()

    return status

def get_pairwise_dist(file_path):
	'''A function that gets the uids and the corresponding scores
	and prints them in tsv.
	'''

	uid_pairs = [] #List with unique ids

	#df.loc[df['pdb'] == '1zmq']
	print('uid1' + '\t' + 'uid2' + '\t' + 'RMSD')
	with open(file_path) as file:
		for line in file:
			if 'Name of Chain' in line:
				line = line.rstrip() #remove \n
				line = line.split("/") #split on /
				uid = line[-1].split(".")[0] #Get uid
				
				uid_pairs.append(uid)

			if 'RMSD=' in line:
				line = line.rstrip() #remove \n
				line = line.split(",") #split on ,
				RMSD = line[1].split(' ')[-1] #Split on space

				print(uid_pairs[0] + '\t' + uid_pairs[1] + '\t' + str(RMSD))
				uid_pairs = [] #reset list of pairs




	return status
#####MAIN#####
args = parser.parse_args()

uid_file = args.uid_file[0]
filter_file = args.filter_file[0]
address = args.address[0]
output_dir = args.output_dir[0]
TMalign = args.TMalign_path[0]

#Get selected uids:
H_group = uid_file.split('/')[-1] #Get H-group (last part of path)

uids = read_newline(uid_file)

#Get pdb ids to filter on
filter_ids = read_newline(filter_file)

downloaded_ids = get_structures(address, uids, filter_ids, H_group, TMalign)

#Print downloaded ids
print(downloaded_ids)
pdb.set_trace()



