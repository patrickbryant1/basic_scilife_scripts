#! /usr/bin/env python3
# -*- coding: utf-8 -*-


import argparse
import sys
import os
import glob
import pdb


#Arguments for argparse module:
parser = argparse.ArgumentParser(description = '''A program that filters uids based on the presence of a parent pdb code in a filter list.''')

parser.add_argument('fasta_file', nargs=1, type= str,
                  default=sys.stdin, help = 'path to .fa file.')
parser.add_argument('filter_file', nargs=1, type= str,
                  default=sys.stdin, help = 'path to file that contains newline separated pdb ids from a pdb search.')
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
			contents.append(line) #Add

	return(contents)



def read_fasta(fasta_file , filter_ids, output_dir):
	'''Read fasta sequences into dict
	'''


	fasta_dict = {} #Store fasta_sequence
	failed_pdb_filter = [] #Store failed


	with open(fasta_file, 'r') as file:
		sequence = ''
		for line in file:
			line = line.rstrip() #remove \n
			if line[0] == '>':
				uid = line.split('|')[-1].split('/')[0]

				if uid[0:4].upper() in filter_ids: #Make check on pdb search
    					continue
				else:
    					failed_pdb_filter.append(uid)

	with open(output_dir+'failed_pdb_filter', 'w') as f:
	               for i in failed_pdb_filter:
	                       f.write(str(i)+'\n')



	return(fasta_dict)

#####MAIN#####
args = parser.parse_args()

fasta_file = args.fasta_file[0]
filter_file = args.filter_file[0]
output_dir = args.output_dir[0]


#Get pdb ids to filter on
filter_ids = read_newline(filter_file)
fasta_dict = read_fasta(fasta_file , filter_ids, output_dir) #Get fasta sequences - filter on filter_ids
pdb.set_trace()
