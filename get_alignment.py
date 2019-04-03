#! /usr/bin/env python3

# -*- coding: utf-8 -*-


import argparse
import sys
import os

import pdb


#Arguments for argparse module:
parser = argparse.ArgumentParser(description = '''A program that takes aligned residue pairs from the structural
								alignment from TMalign and runs tree-puzzle on them.''')
 
parser.add_argument('align_file', nargs=1, type= str,
                  default=sys.stdin, help = 'path to file with TMalign output.')

parser.add_argument('out_dir', nargs=1, type= str,
                  default=sys.stdin, help = 'path to output directory (include / in end).')

#Functions
def get_alignments(file_path, out_dir):
	'''A function that gets the aligned residue pairs from the structural
	   alignment from TMalign and writes them to files in phylip format
	   including copies of each sequence due to the quartet requirement of tree-puzzle.
	'''

	uid_pairs = [] #List with unique ids
	aligned_seqs = [] #List with aligned residue pairs
	fetch_next_3 = False #Keeping track of lines
	file_names = [] #Store file names

	with open(file_path) as file:
		for line in file:
			if 'Name of Chain' in line:
				line = line.rstrip() #remove \n
				line = line.split("/") #split on /
				uid = line[-1].split(".")[0] #Get uid
				
				uid_pairs.append(uid)



			if fetch_next_3 == True: #Fetch next three lines
				#Fetch lines
				if fetched_lines < 3:
					fetched_lines+=1
					line = line.rstrip() #remove \n
					aligned_seqs.append(line) #Append to list
				else:
					file_name = print_alignment(uid_pairs, aligned_seqs, out_dir)
					file_names.append(file_name)
					fetch_next_3 = False
					uid_pairs = [] #reset list of pairs
					aligned_seqs = [] #reset aligned seqs
			
				
			if 'denotes aligned residue pairs' in line:
				fetch_next_3 = True
				fetched_lines = 0
			


	return file_names

def print_alignment(uid_pairs, aligned_seqs, out_dir):
	'''Print alignment for input to RNN
	'''
	
	#Define file name
	file_name = uid_pairs[0] + '_' + uid_pairs[1] + '.aln'
	#Open file and write text to it
	with open(out_dir+file_name, "w") as file:
		file.write(aligned_seqs[0]+'\n')
		file.write(aligned_seqs[2])

	return file_name


#Main program
args = parser.parse_args()

align_file = args.align_file[0]
out_dir = args.out_dir[0]


#Make files with alignments
file_names = get_alignments(align_file, out_dir)
print('Number of extracted per-residue alignments from TMalign: '+str(len(file_names)))

