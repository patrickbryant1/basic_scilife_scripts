#! /usr/bin/env python3
# -*- coding: utf-8 -*-


import argparse
import sys
import os
import matplotlib.pyplot as plt
import numpy
import random

import pdb

#Arguments for argparse module:
parser = argparse.ArgumentParser(description = '''A program that investigates the distribution of
									 each H-group in CATH. And writes all uids for each H-group with at least x entries
									 into files(newline separated) in the output directory''')

parser.add_argument('H_groups_file', nargs=1, type= str,
                  default=sys.stdin, help = 'Path to file with H-groups.')

parser.add_argument('fasta_file', nargs=1, type= str,
                  default=sys.stdin, help = 'Path to fasta file.')

parser.add_argument('outdir_path', nargs=1, type= str,
                  default=sys.stdin, help = 'Path to output directory.')



#Functions
def read_tsv(H_groups_file):
	'''Read ids and H-groups into dict
	'''

	H_groups = {} #Store H-groups and uids

	with open(H_groups_file) as file:
		for line in file:
			line = line.rstrip() #remove \n
			line = line.split(',')
			uid = line[0]
			H_group = line[1:]
			H_group = H_group[0]+'.'+H_group[1]+'.'+H_group[2]+'.'+H_group[3]
			H_groups[uid] = H_group


	return H_groups

def read_fasta(fasta_file):
	'''Read fasta file into dict
	'''

	sequences = {} #Store sequences

	with open(fasta_file) as file:
		for line in file:
			line = line.rstrip() #remove \n
			if line[0] == '>':
				uid = line.split('|')[2].split('/')[0]
			else:
				sequences[uid] = line

	return sequences


def get_groups(H_groups, sequences):
	'''Get H-group for each uid and group sequences accordingly
	'''


	grouped_sequences = {} #Sequences grouped by H-group

	for key in sequences:
		H_group = H_groups[key]
		sequence = sequences[key]

		if H_group not in grouped_sequences.keys(): #If not in new dict - add
			grouped_sequences[H_group] = [key + '/' + sequence]
		else:
			grouped_sequences[H_group].append(key + '/' + sequence) #Otherwise append

	return grouped_sequences

def write_fasta_by_group(grouped_sequences, outdir_path):
	'''Write the selected fasta sequences into directories by H-group
	'''

	for group in grouped_sequences:
		group_dir = outdir_path+group
		os.mkdir(group_dir)
		for fasta in grouped_sequences[group]:
			fasta = fasta.split('/')
			uid = fasta[0]
			sequence = fasta[1]
			with open(group_dir+'/'+uid+'.fa', "w") as file:
				file.write('>'+uid+'\n')
				i = 0 #index
				while i<len(sequence):
					file.write(sequence[i:i+60]+'\n')
					i+=60

	return None


#######################MAIN################################
args = parser.parse_args()

H_groups_file = args.H_groups_file[0]
fasta_file = args.fasta_file[0]
outdir_path = args.outdir_path[0]


H_groups = read_tsv(H_groups_file)
sequences = read_fasta(fasta_file)
grouped_sequences = get_groups(H_groups, sequences)

write_fasta_by_group(grouped_sequences, outdir_path)
