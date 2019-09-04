#! /usr/bin/env python3
# -*- coding: utf-8 -*-


import argparse
import sys
import subprocess
from pathlib import Path
import random
import os
import glob
import pdb

#Custom imports
from conversions import run_hhblits, make_phylip, pdb_to_fasta
from hh_reader import read_result

#Arguments for argparse module:
parser = argparse.ArgumentParser(description = '''A program that converts fasta files to HMMs using hhblits.
						All pairs of x randomly selected entries are then aligned with hhalign.
						''')

parser.add_argument('input_dir', nargs=1, type= str,
                  default=sys.stdin, help = 'path to directory with .fa files.')
parser.add_argument('output_dir', nargs=1, type= str,
                  default=sys.stdin, help = 'output directory.')
parser.add_argument('hhblits', nargs=1, type= str,
                  default=sys.stdin, help = 'path to hhblits.')
parser.add_argument('hhalign', nargs=1, type= str,
                  default=sys.stdin, help = 'path to hhalign.')
parser.add_argument('uniprot', nargs=1, type= str,
                  default=sys.stdin, help = 'path to uniprot database.')
parser.add_argument('get_max', nargs=1, type= int,
                  default=sys.stdin, help = 'Maximum number of structures to use.')
parser.add_argument('address', nargs=1, type= str,
                  default=sys.stdin, help = 'Web adress to download from.') #www.cathdb.info/version/v4_2_0/api/rest/id/

#FUNCTIONS
def read_fasta(input_dir, output_dir):
	'''Read fasta sequences into dict
	'''


	fasta_dict = {} #Store fasta_sequence
	fasta_files = glob.glob(input_dir +'*.fa')


	for file_name in fasta_files:
		with open(file_name) as file:
			sequence = ''
			for line in file:
				line = line.rstrip() #remove \n
				if line[0] == '>':
					uid = line[1:]
				else:
					sequence += line

			fasta_dict[uid] = sequence
	return(fasta_dict)





def loop_through_ids(fasta_dict, uids, H_group,  input_dir, output_dir, hhblits, hhalign, uniprot, get_max, address):
	'''Try to make maximum get_max alignments that fulfill criteria.
	'''

	#Shuffle uids to make sure there is no selective order in comparisons within H-groups
	random.Random(2).shuffle(uids)

	converted_uids = [] #Those uids which .fa have been represented as a HMM
	status = False #Set original status
	identities = {} #Save identities from hhalign

	#Get as many as possible, bu maximum get_max
	get_max = min(get_max, len(uids))
	selected_uids = uids[0:get_max]
	for i in range(len(selected_uids)):
		#Check if .hhm exists
		if not os.path.isfile(output_dir+selected_uids[i]+'.hhm'):
			#Make HMM
			run_hhblits(selected_uids[i], input_dir, hhblits, uniprot)
		#Check if .pdb exists
		if not os.path.isfile(output_dir+selected_uids[i]+'.pdb'):
			#Get pdb file
			subprocess.call(["wget",address+selected_uids[i]+'.pdb'])
	#Align
	(latest_pos, identities) = align(selected_uids, output_dir, H_group, hhalign, identities)

	#Write identities to file
	with open(output_dir+'identities', 'w') as f:
		for key in identities:
			f.write(key+'\t'+str(identities[key])+'\n')

	
	return None

def align(selected_uids, output_dir, H_group, hhalign, identities):
	'''Run hhalign on file pair and extract sequence identity and
	% aligned of shortest sequence.
	'''



	end = len(selected_uids)

	parsed_output = {} #Save parsed output from hhalign
	for i in range(0, end):
		structure_i =output_dir+selected_uids[i]+'.hhm' #Get structure i

		for j in range(i+1, end):
			structure_j =output_dir+selected_uids[j]+'.hhm' #Get structure j
			#Run hhalign
			outp = subprocess.check_output(hhalign + ' -i '+ structure_i + ' -t ' + structure_j + ' -o ' + uids[i]+'_'+uids[j]+'.hhr' + ' -norealign', shell=True)
			result = read_result(output_dir+uids[i]+'_'+uids[j]+'.hhr')
			chain_lens = [result[0].query_length, result[0].template_length]
			aligned_len = result[0].aligned_cols
			identity = result[0].identity
			query_aln = result[0].query_ali
			template_aln = result[0].template_ali
			start_pos = result[0].start
			end_pos = result[0].end
			e_value = result[0].evalue
			probability = result[0].probability
			#Save identities to see distributions
			key = selected_uids[i]+'_'+selected_uids[j]
			if key not in identities.keys():
				identities[key] = identity

			#Add to write			
			parsed_output[str(selected_uids[i]+'_'+selected_uids[j])] = [query_aln, template_aln, chain_lens, aligned_len, identity, start_pos, end_pos, e_value, probability] #Add info to parsed output


	write_to_file(output_dir, H_group, parsed_output)

	return(i, identities)


def write_to_file(output_dir, H_group, parsed_output):
	'''Write all extracted information from hhalign to files
	that will be used downstream.
	'''

	for key in parsed_output:
		#Get uids and saved aln info
		uids = key.split('_')
		query_aln, template_aln, chain_lens, aligned_len, identity, start_pos, end_pos, e_value, probability = parsed_output[key]
		#Write alignment and info to file
		with open(output_dir+key+'.aln', 'w') as f:
			#f.write('#'+'query:' + 'l=' + str(chain_lens[0]) + ' s=' + str(start_pos[0]) + ' e=' + str(end_pos[0]) + '|template: ' + 'l=' + str(chain_lens[1]) + ' s=' + str(start_pos[1]) + ' e=' + str(end_pos[1]) +  '|aligned_len: ' + str(aligned_len) + '|Identity: ' + str(identity) + '\n')
			f.write('>'+uids[0]+'|l='+str(chain_lens[0]) + '|s=' + str(start_pos[0]) + '|e=' + str(end_pos[0])+'|aligned_len=' + str(aligned_len) + '|Identity=' + str(identity) + 
			'|e-value='+ str(e_value) + '|probability=' + str(probability) + '\n')
			f.write(query_aln+'\n') #write sequences
			f.write('>'+uids[1]+'|l=' + str(chain_lens[1]) + '|s=' + str(start_pos[1]) + '|e=' + str(end_pos[1])+'\n')
			f.write(template_aln)

		#Write .phy file of alignment
		make_phylip(uids, query_aln, template_aln, output_dir)

	return None

#####MAIN#####
args = parser.parse_args()

input_dir = args.input_dir[0]
output_dir = args.output_dir[0]
hhblits = args.hhblits[0]
hhalign = args.hhalign[0]
uniprot = args.uniprot[0]
get_max = args.get_max[0]
address = args.address[0]
H_group = input_dir.split('/')[-1] #Get H-group (last part of path)

#Get pdb ids to filter on
fasta_dict = read_fasta(input_dir, output_dir) #Get fasta sequences - filter on filter_ids
uids = [*fasta_dict.keys()] #Get uids
loop_through_ids(fasta_dict, uids, H_group, input_dir,  output_dir, hhblits, hhalign, uniprot, get_max, address)

