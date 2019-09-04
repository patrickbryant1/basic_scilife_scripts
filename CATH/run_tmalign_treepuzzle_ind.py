#! /usr/bin/env python3
# -*- coding: utf-8 -*-


import argparse
import sys
import os
import glob
import subprocess
from conversions import make_phylip 
from match_pdb import seq_to_pdb
import pdb



#Arguments for argparse module:
parser = argparse.ArgumentParser(description = '''A program that runs TMalign and tree-puzzle and
						receives the resulting output.''')

parser.add_argument('indir', nargs=1, type= str, default=sys.stdin, help = 'Path to input directory.')
parser.add_argument('outdir', nargs=1, type= str, default=sys.stdin, help = 'Path to output directory.')
parser.add_argument('fastadir', nargs=1, type= str, default=sys.stdin, help = 'path to directory with .fa files.')
parser.add_argument('hgroup', nargs=1, type= str, default=sys.stdin, help = 'H-group.')
parser.add_argument('puzzle', nargs=1, type= str, default=sys.stdin, help = 'Path to tree-puzzle.')
parser.add_argument('TMalign', nargs=1, type= str, default=sys.stdin, help = 'Path to TMalign.')


#FUNCTIONS
def read_fasta(filepath):
	'''Reads aligned sequences in fasta format
	'''
	
	fasta_dict = {} #Store fasta_sequence
	with open(filepath, 'r') as file:
		sequence = ''
		fetched = False
		for line in file:
			line = line.rstrip()
			if line[0] == '>':
				if fetched == True:
					fasta_dict[uid] = sequence
					sequence = '' #Reset sequence
					fetched = False
				uid = line[1:8]
			else:
				sequence += line
				fetched = True
			
	#Add last sequence
	fasta_dict[uid] = sequence
	return fasta_dict

def run_puzzle(outdir, puzzle):
	'''Run tree-puzzle and retrieve output
	'''
	for name in glob.glob(outdir+"*.phy"): #Use all .phy files
		uid_pairs = name.split('/')[-1].split('.')[0].split('_')
		try:
			p = subprocess.Popen([puzzle, name], stdin=subprocess.PIPE)
			p.communicate(b'y\nn\n')[0]
		except:
			raise IOError(name)


	return None

def run_TMalign(indir, outdir, fastadir, TMalign):
	'''Run TMalign on .pdb files
	'''
	
	measures = {} #Save RMSD to add with MLAA distance from tree-puzzle
	names = glob.glob(indir+"*.aln") #Use all .aln files
	status = True #See if H-group has enough entries fulfilling criteria
	n = 1 #at least n structures compared
	if len(names) < (n):
		status = False
	if status == True:
		while names:#While names not empty
			aln_i = names[0] #Get structure i
			uids = aln_i.split('/')[-1].split('.')[0].split('_')
			uid1 = uids[0]
			uid2 = uids[1]
			names.pop(0)
			#Run TMalign and extract scores
			str1 = indir+uid1+'.pdb'
			str2 = indir+uid2+'.pdb'

			tmalign_out = subprocess.check_output([TMalign, str1 , str2]) #Performs optimal structural alignment 
			(tm_aligned_len, rmsd, tmscores, tm_identity, chain_lens, tm_sequences)= parse_tm(tmalign_out)	
			measures[uid1+'_'+uid2] = [rmsd, tmscores[0], tmscores[1]]
			#Write .phy file of alignment
			make_phylip(uids, tm_sequences[0], tm_sequences[1], outdir)
			#Get original fasta sequences
			org1 = read_fasta(fastadir+uid1+'.fa')
			org2 = read_fasta(fastadir+uid2+'.fa')
			#Write new .pdb files matching alignment to be used for lddt
			seq_to_pdb(uids, tm_sequences[0], tm_sequences[1], org1[uid1], org2[uid2], outdir)
			#Write the alignment
			with open(outdir+uid1+'_'+uid2+'.aln', 'w') as f:
				f.write('>'+uids[0]+'|l='+str(chain_lens[0]) + '|aligned_len=' + str(tm_aligned_len) + '|Identity=' + str(tm_identity)+'\n')
				f.write(tm_sequences[0]+'\n') #write sequences
				f.write('>'+uids[1]+'|l=' + str(chain_lens[1])+'\n')
				f.write(tm_sequences[1])

	return measures, status

def parse_tm(tmalign_out):
	'''A function that gets the uids and the corresponding scores
	and prints them in tsv.
	'''
	
	tmalign_out = tmalign_out.decode("utf-8")
	tmalign_out = tmalign_out.split('\n')
	tmscores = [] #Save TMscores
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
		if 'TM-score=' in tmalign_out[i]:
			tmscores.append(tmalign_out[i].split('(')[0].split('=')[1].strip())

	#Get per residue sequence alignments from structural alignment
	sequences = [tmalign_out[-5], tmalign_out[-3]]

	chain_lens = [int(len_1), int(len_2)]
	
			
	return(aligned_len, rmsd, tmscores, identity, chain_lens, sequences)


def parse_puzzle(measures, indir):
	'''Parse output from tree-puzzle and write to dict
	'''
	keys = [*measures] #Make list of keys in dict
	for key in keys:
		uids = key.split('_')
		rmsd, tmscore1, tmscore2 = measures[key] #Get rmsd
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
				measures[key] = [rmsd, tmscore1, tmscore2, seq_dist] 
				break
		dist_file.close()

	return measures


def print_tsv(measures, hgroup, outdir):
	'''Print measures in tsv to file
	'''
	with open(outdir+hgroup+'_str.tsv', 'w') as file:
		file.write('uid1\tuid2\tMLAAdist\tRMSD\tTMscore_high\tTMscore_low\n')
		for key in measures:
			uids = key.split('_')
			rmsd, tmscore1, tmscore2, seq_dist = measures[key]
			high_score = max(float(tmscore1), float(tmscore2))
			low_score = min(float(tmscore1), float(tmscore2)) 
			file.write(uids[0]+'\t'+uids[1]+'\t'+seq_dist+'\t'+rmsd+'\t'+str(high_score)+'\t'+str(low_score)+'\n')

	return None

#####MAIN#####
args = parser.parse_args()

indir = args.indir[0]
outdir = args.outdir[0]
fastadir = args.fastadir[0]
hgroup = args.hgroup[0]
puzzle = args.puzzle[0]
TMalign = args.TMalign[0]

measures, status = run_TMalign(indir, outdir, fastadir, TMalign)
run_puzzle(outdir, puzzle)

if status == True: #Only if H-groups fulfills criteria
	measures = parse_puzzle(measures, outdir)
	print_tsv(measures, hgroup, outdir)
