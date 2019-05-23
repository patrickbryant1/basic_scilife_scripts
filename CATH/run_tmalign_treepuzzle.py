#! /usr/bin/env python3
# -*- coding: utf-8 -*-


import argparse
import sys
import os
import glob
import subprocess
import pdb



#Arguments for argparse module:
parser = argparse.ArgumentParser(description = '''A program that runs TMalign and tree-puzzle and
						receives the resulting output.''')

parser.add_argument('indir', nargs=1, type= str, default=sys.stdin, help = 'Path to input directory.')
parser.add_argument('puzzle', nargs=1, type= str, default=sys.stdin, help = 'Path to tree-puzzle.')
parser.add_argument('TMalign', nargs=1, type= str, default=sys.stdin, help = 'Path to TMalign.')


#FUNCTIONS


def run_puzzle(indir, puzzle):
	'''Run tree-puzzle and retrieve output
	'''
	for name in glob.glob(indir+"*.phy"): #Use all .phy files
		uid_pairs = name.split('/')[-1].split('.')[0].split('_')
		try:
			p = subprocess.Popen([puzzle, name], stdin=subprocess.PIPE)
			p.communicate(b'y\nn\n')[0]
		except:
			raise IOError(name)


	return None

def run_TMalign(indir, TMalign):
	'''Run TMalign on extracted CAs from hhalign alignments
	'''
	
	measures = {} #Save RMSD to add with MLAA distance from tree-puzzle
	names = glob.glob(indir+"*_aln.pdb") #Use all _aln.pdb files
	
	for i in range(0, len(names)):
		structure_i = names[i] #Get structure i
		uid1 = structure_i.split('/')[-1].split('_')[0]
		for j in range(i+1, len(names)):
			structure_j =names[j] #Get structure j
			uid2 = structure_j.split('/')[-1].split('_')[0]
			tmalign_out = subprocess.check_output([TMalign, structure_i , structure_j , '-a'])
			(tm_aligned_len, rmsd, tm_identity, chain_lens, tm_sequences)= parse_tm(tmalign_out)	
			measures[uid1+'_'+uid2] = [rmsd]
	return measures

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


def parse_puzzle(measures, indir):
	'''Parse output from tree-puzzle and write to dict
	'''
	for key in measures:
		name =indir + key + '.phy.dist'
		dist_file = open(name, 'r') 
		for line in dist_file:
			line = line.rstrip()
			line = line.split(" ") #split on double space
			line = list(filter(None, line)) #Filter away empty strings

			if len(line)>2:
				seq_dist = line[-1] #Get ML evolutionary distance between sequences
				pdb.set_trace()
				measures[key] = measures[key].append(seq_dist) 
				break


	return measures


#####MAIN#####
args = parser.parse_args()

indir = args.indir[0]
puzzle = args.puzzle[0]
TMalign = args.TMalign[0]

run_puzzle(indir, puzzle)
measures = run_TMalign(indir, TMalign)
measures = parse_puzzle(measures, indir)
pdb.set_trace()
