#! /usr/bin/env python3

# -*- coding: utf-8 -*-


import argparse
import sys
import os
import subprocess
import glob 
import pdb


#Arguments for argparse module:
parser = argparse.ArgumentParser(description = '''A program that takes aligned residue pairs from the structural
								alignment from TMalign and runs tree-puzzle on them.''')
 
parser.add_argument('align_dir', nargs=1, type= str,
                  default=sys.stdin, help = 'Path to directory with TMalign alignments.')

#Functions
def get_alignments(align_dir):
	'''A function that gets the aligned residue pairs from the structural
	   alignment from TMalign and writes them to files in phylip format
	   including copies of each sequence due to the quartet requirement of tree-puzzle.
	'''

	
	file_names = [] #Store file names
	
	for infile in glob.glob(align_dir+"*.aln"):
		
		with open(infile) as file:
			infile = infile.split('/')[-1] #Last part of path
			uid_pairs =  infile.split('.')[0].split('_')
			aligned_seqs = [] #Store aligned sequences
			#Get sequence alignments			
			for line in file:
				aligned_seqs.append(line.rstrip())
			#Format sequence alignment into phylip to run tree-puzzle
			file_name = make_phylip(align_dir, uid_pairs, aligned_seqs)
			file_names.append(file_name)
						

	
	return file_names

def make_phylip(align_dir, uid_pairs, aligned_seqs):
	'''Print phylip format for tree-puzzle calculations
	'''
	#Create text in phylip format
	text = (' 4  ' + str(len(aligned_seqs[0])) + '\n'
			+ uid_pairs[0] + '00|' + aligned_seqs[0] + '\n'
			+ 'copy11111' + '|' + aligned_seqs[0] + '\n'
			+ uid_pairs[1] + '00|' + aligned_seqs[1] + '\n'
			+ 'copy22222' + '|' + aligned_seqs[1] + '\n')
	

	#Define file name
	file_name = align_dir+uid_pairs[0] + '_' + uid_pairs[1] + '.phy'
	#Open file and write text to it
	with open(file_name, "w") as file:
		file.write(text)

	return file_name


#Main program
args = parser.parse_args()

align_dir = args.align_dir[0]



#Make .phy files with alignments
file_names = get_alignments(align_dir)

print('Number of files to tree-puzzle: '+str(len(file_names)))
#Run tree-puzzle on the files
for name in file_names:
	message_1 = "/home/p/pbryant/pfs/tree-puzzle-5.3.rc16-linux/src/puzzle"
	try:
		p = subprocess.Popen([message_1, name], stdin=subprocess.PIPE)
		p.communicate(b'y\nn\n')[0]
		
	except:
		raise IOError(name)
