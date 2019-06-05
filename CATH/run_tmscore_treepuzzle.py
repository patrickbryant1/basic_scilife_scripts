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
parser.add_argument('hgroup', nargs=1, type= str, default=sys.stdin, help = 'H-group.')
parser.add_argument('puzzle', nargs=1, type= str, default=sys.stdin, help = 'Path to tree-puzzle.')
parser.add_argument('TMscore', nargs=1, type= str, default=sys.stdin, help = 'Path to TMscore.')


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

def run_TMscore(indir, TMscore):
	'''Run TMalign on extracted CAs from hhalign alignments
	'''
	
	measures = {} #Save RMSD to add with MLAA distance from tree-puzzle
	names = glob.glob(indir+"*_aln.pdb") #Use all _aln.pdb files
	done_uids = [] #Keep trackof the uids that are completed
	status = True #See if H-group has enough entries fulfilling criteria
	n = 5 #at least n structures compared
	if len(names) < (n*2):
		status = False
	if status == True:
		while names:#While names not empty
			structure_i = names[0] #Get structure i
			uid1 = structure_i.split('/')[-1].split('_')[0]
			uid2 = structure_i.split('/')[-1].split('_')[2]
			names.pop(0)
			for j in range(0, len(names)):
				structure_j =names[j] #Get structure j
				if uid1 in structure_j and uid2 in structure_j: #If uid1 and uid2 is part of the file,
										#it is the right pdb representation of the alignment
					tscore_out = subprocess.check_output([TMscore, structure_i , structure_j , '-a'])
					(tm_aligned_len, rmsd, tm_identity, chain_lens, tm_sequences)= parse_tm(tmscore_out)	
					measures[uid1+'_'+uid2] = rmsd
					break #Break, since match found
					names.pop(j)

	return measures, status

def parse_tm(tmscore_out):
	'''A function that gets the uids and the corresponding scores
	and prints them in tsv.
	'''
	
	tmscore_out = tmscore_out.decode("utf-8")
	tmscore_out = tmscore_out.split('\n')
	
	for i in range(0, len(tmscore_out)): #Step through all items in list
		if 'TM-score    =' in tmscore_out[i]:
			row = tmscore_out[i].split('=')
                        score = row[1].split('(')[0].rstrip() 
		if 'Superposition in the TM-score:' in tmalign_out[i]:
			row = tmscore_out[i].split('=')
			rmsd = row[-1].rstrip()	
				

	pdb.set_trace()
	
			
	return(score, rmsd)


def parse_puzzle(measures, indir):
	'''Parse output from tree-puzzle and write to dict
	'''
	keys = [*measures] #Make list of keys in dict
	for key in keys:
		uids = key.split('_')
		rmsd = measures[key] #Get rmsd
		try:
			dist_file = open(indir + key + '.phy.dist')
		except:
			uids = key.split('_')
			dist_file = open(indir + uids[1] + '_' + uids[0] + '.phy.dist')
			measures.pop(key)
			#change key to match other file names
			key = uids[1] + '_' + uids[0]
		for line in dist_file:
			line = line.rstrip()
			line = line.split(" ") #split on double space
			line = list(filter(None, line)) #Filter away empty strings

			if len(line)>2:
				seq_dist = line[-1] #Get ML evolutionary distance between sequences
				measures[key] = [rmsd, seq_dist] 
				break
		dist_file.close()

	return measures


def print_tsv(measures, hgroup):
	'''Print measures in tsv to file
	'''
	with open(hgroup+'.tsv', 'w') as file:
		file.write('uid1\tuid2\tMLAAdist\tRMSD\n')
		for key in measures:
			uids = key.split('_')
			info = measures[key]
			rmsd, seq_dist = info[0], info[1]
			file.write(uids[0]+'\t'+uids[1]+'\t'+seq_dist+'\t'+rmsd+'\n')

	return None

#####MAIN#####
args = parser.parse_args()

indir = args.indir[0]
hgroup = args.hgroup[0]
puzzle = args.puzzle[0]
TMscore = args.TMscore[0]

run_puzzle(indir, puzzle)
(measures, status) = run_TMscore(indir, TMscore)
if status == True: #Only if H-groups fulfills criteria
	measures = parse_puzzle(measures, indir)
	print_tsv(measures, hgroup)
