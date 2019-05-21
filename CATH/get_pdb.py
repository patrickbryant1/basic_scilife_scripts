#! /usr/bin/env python2
# -*- coding: utf-8 -*-


import argparse
import sys
import subprocess
import pexpect
import random
import os
import pdb

#Custom imports
from hh_suite import pdb_to_fasta, run_hhblits
from hh_reader import read_result

#Arguments for argparse module:
parser = argparse.ArgumentParser(description = '''A program that downloads pdb structures based on CATH uids (domain ids)
						from H-groups that have at least x entries.
						It then converts the files to fasta and to a HMM using hhblits.
						All pairs of x-y randomly selected entries are then aligned with hhalign. 
						If less than 75 % of the shortest sequence has been aligned, the second uid
						is dropped and a new domain structure downloaded, if available.
						If the aligned sequences are over 90 % identical, the second uid
                                                is dropped and a new domain structure downloaded, if available. 
						Otherwise the whole H-group is dropped''')
 
parser.add_argument('uid_file', nargs=1, type= str,
                  default=sys.stdin, help = 'path to CATH ids file.')
parser.add_argument('filter_file', nargs=1, type= str,
                  default=sys.stdin, help = 'path to file that contains newline separated pdb ids from a pdb search.')
parser.add_argument('address', nargs=1, type= str,
                  default=sys.stdin, help = 'Web adress to download from.') #www.cathdb.info/version/v4_2_0/api/rest/id/
parser.add_argument('output_dir', nargs=1, type= str,
                  default=sys.stdin, help = 'output directory.')
parser.add_argument('hhblits', nargs=1, type= str,
                  default=sys.stdin, help = 'path to hhblits.')
parser.add_argument('hhalign', nargs=1, type= str,
                  default=sys.stdin, help = 'path to hhalign.')
parser.add_argument('uniprot', nargs=1, type= str,
                  default=sys.stdin, help = 'path to uniprot database.')

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



def get_structures(address, uids, filter_ids, H_group, output_dir, hhblits, hhalign, uniprot):
	'''Download .pdb structure, run TMalign and check sequence identity
	'''

	downloaded_uids = [] #Keep track of ids that have been downloaded	

	#status = False #Set original status
	#Shuffle uids to make sure there is no selective order in comparisons within H-groups
	random.Random(2).shuffle(uids)
	
	#Check how many uids to get
	if len(uids) > 5:
		get_n = 5
	else:
		get_n = len(uids)

	#Go through uids and try to find get_n uids that match criteria
	while get_n >= 2:
		(status, downloaded_uids, selected_uids) = loop_through_ids(address, uids, filter_ids, H_group, output_dir, get_n, downloaded_uids, hhblits, hhalign, uniprot)
		if status == True:
			#Remove the uids not to be used for later steps
			for duid in downloaded_uids:
				if duid not in selected_uids:
					os.remove(output_dir + duid +'.pdb') #Remove failed .pdb file
			break
		else:
			get_n -= 1
		
	#If you could not get at least 2 uids that fulfill criteria
	if status == False:
		print('The H-group ' + H_group + ' does not fulfill the criteria.')

	print(downloaded_uids)
        
	return None

def loop_through_ids(address, uids, filter_ids, H_group, output_dir, get_n, downloaded_uids, hhblits, hhalign, uniprot):
	'''Loop through uids and try to make get_n alignments
	that fulfill criteria.
	'''
	selected_uids = [] #Selected for hhsuite
	status = False #Set original status

	for i in range(0, len(uids)):
		if len(selected_uids) == get_n:
		#Make alignment of all of these

			(status, latest_pos) = align(selected_uids, output_dir, H_group, hhalign)
			#If one fails, pop this and continue
			if status == False: #If all have not passed
				selected_uids.pop(latest_pos) #Pop this failed uid
			else:
				print(H_group + 'aligned!')
				break


		if uids[i][0:4].upper() in filter_ids: #Make check on pdb search
			selected_uids.append(uids[i])
			if uids[i] not in downloaded_uids:
				downloaded_uids.append(uids[i])
				#Get pdb file
				subprocess.call(["wget",address+uids[i]+'.pdb'])
				#Make fasta
				pdb_to_fasta(uids[i], output_dir)
				#Make HMM
				run_hhblits(uids[i], output_dir, hhblits, uniprot)
				
		else:
			print(uids[i] + ' did not pass filter')

	return(status, downloaded_uids, selected_uids)

def align(selected_uids, output_dir, H_group, hhalign):
	'''Run hhalign on file pair and extract sequence identity and
	% aligned of shortest sequence.
	Remove file2 if identity is 90 % or above or if less than 75 % has been aligned.
	'''


    
	count = 0 #Keep track of number of alignments made
	end = len(selected_uids)
	status = True #Keep track on if sequences are too similar
		      #If True - no problem

	

	parsed_output = {} #Save parsed output from hhalign

	for i in range(0, end):
		structure_i =output_dir+selected_uids[i]+'.hhm' #Get structure i
		
		for j in range(i+1, end):
			structure_j =output_dir+selected_uids[j]+'.hhm' #Get structure j
			#Run hhalign
			pexpect.run(hhalign + ' -i '+ structure_i + ' -t ' + structure_j + ' -o ' + uids[i]+'_'+uids[j]+'.hhr' + ' -glob')
			result = read_result(output_dir+uids[i]+'_'+uids[j]+'.hhr')
			chain_lens = [result[0].query_length, result[0].template_length]
			aligned_len = result[0].aligned_cols
			identity = result[0].identity
			query_aln = result[0].query_ali
			template_aln = result[0].template_ali
			
			 
			if (aligned_len < (0.75*min(chain_lens))) or (identity >= 0.90): #aligned lenght and sequence identity thresholds 
				print(selected_uids[i], selected_uids[j])
				print(aligned_len, shortest_seq, identity)
				status = False
				break #Break out, since too similar seqs
			count+=1    
			parsed_output[str(selected_uids[i]+'_'+selected_uids[j])] = [query_aln, template_aln, chain_lens, aligned_len, identity] #Add info to parsed output
			
		if status == False:
			break #Break out, since too similar seqs
			
				 
	if status == True: # save all info to files
				
		write_to_file(output_dir, H_group, parsed_output)
	
	
	return(status, i)


def write_to_file(output_dir, H_group, parsed_output):
	'''Write all extracted information from hhalign to files 
	that will be used downstream
	'''

	for key in parsed_output:
		with open(output_dir+key+'.aln', 'w') as f:
			f.write('#'+'q_length: ' + str(parsed_output[key][2][0]) + '|t_length: ' + str(parsed_output[key][2][1]) + '|aligned_len: ' + str(parsed_output[key][3]) + '|Identity: ' + str(parsed_output[key][4]) + '\n') 
			f.write(parsed_output[key][0]+'\n') #write sequences
                        f.write(parsed_output[key][1])




	return None

#####MAIN#####
args = parser.parse_args()

uid_file = args.uid_file[0]
filter_file = args.filter_file[0]
address = args.address[0]
output_dir = args.output_dir[0]
hhblits = args.hhblits[0]
hhalign = args.hhalign[0]
uniprot = args.uniprot[0]

#Get selected uids:
H_group = uid_file.split('/')[-1] #Get H-group (last part of path)

uids = read_newline(uid_file)

#Get pdb ids to filter on
filter_ids = read_newline(filter_file)

get_structures(address, uids, filter_ids, H_group, output_dir, hhblits, hhalign, uniprot)




