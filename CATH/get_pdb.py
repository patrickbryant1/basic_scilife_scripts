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



def get_structures(address, uids, filter_ids, H_group, TMalign, output_dir):
	'''Download .pdb structure, run TMalign and check sequence identity
	'''

	downloaded_uids = [] #Keep track of ids that have been downloaded	
	selected_uids = [] #Selected for TMalign
	for i in range(0, len(uids)):
		if len(selected_uids) == 10:
			print(selected_uids)
			#Make alignment of all of these
			#pdb.set_trace()
			(status, last_pos) = align(selected_uids, TMalign, output_dir, H_group)
			#If one fails, pop this and continue
			if status == False: #If all have not passed
				selected_uids.pop(last_pos) #Pop this failed uid
			else:
				print('Done')
				break				
			#Save the passed uids to a special list

		if uids[i][0:4].upper() in filter_ids: #Make check on pdb search
			
			downloaded_uids.append(uids[i])
			selected_uids.append(uids[i])
			print(uids[i])
			subprocess.call(["wget",address+uids[i]+'.pdb'])
            
		

	print(downloaded_uids)
        
	return None


def align(selected_uids, TMalign, output_dir,  H_group):
	'''Run TMalign on file pair and extract sequence identity,
	remove file2 if 90 % or above
	'''


    
	count = 0 #Keep track of number of alignments made
	end = len(selected_uids)
	status = True #Keep track on if sequences are too similar
		      #If True - no problem

	

	parsed_output = {} #Save parsed output from TMalign

	for i in range(0, end):
		structure_i =output_dir+selected_uids[i]+'.pdb' #Get structure i
		
		for j in range(i+1, end):
			structure_j =output_dir+selected_uids[j]+'.pdb' #Get structure j
			tmalign_out = subprocess.check_output([TMalign, structure_i , structure_j , '-a'])
			(aligned_len, rmsd, identity, chain_lens, sequences)= parse_tm(tmalign_out)
			if int(aligned_len) > (0.9*min(chain_lens)) and float(identity) > 0.90: #seq identity threshold. Considers the length of the alignment 
				print(selected_uids[i], selected_uids[j])
				print(aligned_len, min(chain_lens), identity)
				status = False
				break #Break out, since too similar seqs
			count+=1    
			parsed_output[str(selected_uids[i]+'_'+selected_uids[j])] = [sequences, rmsd] #Add info to parsed output
			
		if status == False:
			break #Break out, since too similar seqs
			
				 
	if status == True: # save all info to files
		pdb.set_trace()
		write_to_file(output_dir, H_group, parsed_output)
	
	
	return(status, i)

def parse_tm(tmalign_out):
	'''A function that gets the uids and the corresponding scores
	and prints them in tsv.
	'''
	
	tmalign_out = tmalign_out.decode("utf-8")
	tmalign_out = tmalign_out.split('\n')
	
	for i in range(0, len(tmalign_out)): #Step through all items in list
			
		if 'Aligned length' and 'RMSD' and 'Seq_ID' in tmalign_out[i]:
			row = tmalign_out[i].split(',')
			aligned_len = row[0].split('=')[1].lstrip()
			rmsd = row[1].split('=')[1].lstrip()
			identity = row[2].split('=')[2].lstrip() 
		
		if 'Length of Chain_1:' in tmalign_out[i]:
			len_1 = tmalign_out[i].split(':')[1].split()[0]
				
		if 'Length of Chain_2:' in tmalign_out[i]:
                        len_2 = tmalign_out[i].split(':')[1].split()[0]


	#Get per residue sequence alignments from structural alignment
	sequences = [tmalign_out[-5], tmalign_out[-3]]

	chain_lens = [int(len_1), int(len_2)]
	
			
	return(aligned_len, rmsd, identity, chain_lens, sequences)

def write_to_file(output_dir, H_group, parsed_output):
	'''Write all extracted information from TMalign to files 
	that will be used downstream
	'''

	with open(output_dir+H_group+'_rmsd.tsv', 'w') as file_1:
		file_1.write('uid1' + '\t' + 'uid2' + '\t' + 'RMSD') #Write headers
		for key in parsed_output:
			uid_1 = key.split('_')[0]
			uid_2 = key.split('_')[1]
			
			file_1.write(uid_1 + '\t' + uid_2 + '\t' + parsed_output[key][1] + '\n') #Write uids and corresponding rmsd
			with open(output_dir+key+'.aln', 'w') as file_2:
				file_2.write(parsed_output[key][0][0]+'\n') #write sequences
				file_2.write(parsed_output[key][0][1])



	return None

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

downloaded_ids = get_structures(address, uids, filter_ids, H_group, TMalign, output_dir)

#Print downloaded ids
print(downloaded_ids)
pdb.set_trace()



