#! /usr/bin/env python3
# -*- coding: utf-8 -*-


import argparse
import sys
import os
import glob
import subprocess
import pdb

from match_pdb import seq_to_pdb
#Arguments for argparse module:
parser = argparse.ArgumentParser(description = '''A program that runs TMalign and tree-puzzle and
						receives the resulting output.''')

parser.add_argument('indir', nargs=1, type= str, default=sys.stdin, help = 'Path to input directory.')
parser.add_argument('outdir', nargs=1, type= str, default=sys.stdin, help = 'Path to output directory.')
parser.add_argument('fastadir', nargs=1, type= str, default=sys.stdin, help = 'path to directory with .fa files.')
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

def read_fasta(aln_file):
	'''Reads aligned sequences in fasta format
	'''
	
	fasta_dict = {} #Store fasta_sequence
	with open(aln_file) as file:
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

def run_TMscore(indir, fastadir, TMscore):
	'''Run TMalign on extracted CAs from hhalign alignments
	'''
	
	measures = {} #Save RMSD to add with MLAA distance from tree-puzzle
	names = glob.glob(indir+"*.aln") #Use all _aln files
	done_uids = [] #Keep trackof the uids that are completed
	status = True #See if H-group has enough entries fulfilling criteria
	n = 1 #at least n structures compared
	if len(names) < n:
		status = False
	if status == True:
		while names:#While names not empty
			aln_i = names[0] #Get structure i
			uids = aln_i.split('/')[-1].split('.')[0].split('_') 
			uid1 = uids[0]
			uid2 = uids[1]
			
			sequences = read_fasta(aln_i)	
			#Get original fasta sequences
			org1 = read_fasta(fastadir+uid1+'.fa')
			org2 = read_fasta(fastadir+uid2+'.fa')
			#Write new .pdb files matching alignment
			seq_to_pdb(uids, sequences[uid1], sequences[uid2], org1[uid1], org2[uid2] , indir)
			structure_1 = indir+uid1+'_to_'+uid2+'_aln.pdb'
			structure_2 = indir+uid2+'_to_'+uid1+'_aln.pdb'
			try:
				tmscore_out = subprocess.check_output([TMscore, structure_1 , structure_2])
			except:
				names.pop(0) #remove since done
				continue #The files could not be aligned 

			(rmsd, tmscore, gdt_ts, gdt_ha)= parse_tm(tmscore_out)
			if not rmsd:
				print('No common residues btw ' + structure_1 + ' and ' + structure_2 + '\n')	
			else:
				measures[uid1+'_'+uid2] = [rmsd, tmscore, gdt_ts, gdt_ha]
		
			names.pop(0) #remove since done

	return measures, status

def parse_tm(tmscore_out):
	'''A function that gets the uids and the corresponding scores
	and prints them in tsv.
	'''
	
	tmscore_out = tmscore_out.decode("utf-8")
	tmscore_out = tmscore_out.split('\n')
	rmsd = ''
	tmscore = ''
	gdt_ts = ''
	gdt_ha = ''
	for i in range(0, len(tmscore_out)): #Step through all items in list
		if 'TM-score    =' in tmscore_out[i]:
			row = tmscore_out[i].split('=')
			tmscore = row[1].split('(')[0].strip()
		if 'Superposition in the TM-score:' in tmscore_out[i]:
			row = tmscore_out[i].split('=')
			rmsd = row[-1].strip()	
		if 'GDT-TS-score' in tmscore_out[i]:
			gdt_ts = tmscore_out[i].split('%')[0].split('=')[1].strip()		
		if 'GDT-HA-score' in tmscore_out[i]:
			gdt_ha = tmscore_out[i].split('%')[0].split('=')[1].strip()
		if 'There are no common residues in the input structures' in tmscore_out[i]:
			break
						
	return(rmsd, tmscore, gdt_ts, gdt_ha)


def parse_puzzle(measures, indir):
	'''Parse output from tree-puzzle and write to dict
	'''
	keys = [*measures] #Make list of keys in dict
	for key in keys:
		uids = key.split('_')
		rmsd, tmscore, gdt_ts, gdt_ha =  measures[key]	
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
				measures[key] = [seq_dist, rmsd, tmscore, gdt_ts, gdt_ha]
				break
		dist_file.close()
	return measures


def print_tsv(measures, hgroup, outdir):
	'''Print measures in tsv to file
	'''
	with open(outdir+hgroup+'_seq.tsv', 'w') as file:
		file.write('uid1\tuid2\tMLAAdist\tRMSD\tTMscore\tGDT-TS\tGDT-HA\n')
		for key in measures:
			uids = key.split('_')
			info = measures[key]
			seq_dist, rmsd, tmscore, gdt_ts, gdt_ha = info
			file.write(uids[0]+'\t'+uids[1]+'\t'+seq_dist+'\t'+rmsd+'\t'+tmscore+'\t'+gdt_ts+'\t'+gdt_ha+'\n')

	return None

#####MAIN#####
args = parser.parse_args()

indir = args.indir[0]
outdir = args.outdir[0]
fastadir = args.fastadir[0]
hgroup = args.hgroup[0]
puzzle = args.puzzle[0]
TMscore = args.TMscore[0]

run_puzzle(indir, puzzle)
(measures, status) = run_TMscore(indir, fastadir, TMscore)
if status == True: #Only if H-groups fulfills criteria
	measures = parse_puzzle(measures, indir)
	print_tsv(measures, hgroup, outdir)
