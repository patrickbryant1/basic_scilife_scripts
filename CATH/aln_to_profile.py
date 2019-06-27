#! /usr/bin/env python3
# -*- coding: utf-8 -*-


import argparse
import sys
import glob
import subprocess
import pandas as pd
#Custom import
from parse_hmm import read_hmm
#Debugging
import pdb


#Arguments for argparse module:
parser = argparse.ArgumentParser(description = '''A program that matches a fasta alignment to HMMs and 
						returns parsed emission probabilities written to a .csv file.''')

parser.add_argument('indir', nargs=1, type= str,
default=sys.stdin, help = 'path to directory with files.')


#FUNCTIONS
def match_hh(indir):
	'''Match alignments to original HMMs and parse output to csv'''
	
	#Save matched hmms
	matched_hmms =  pd.DataFrame(columns = ['uid1','uid2','hmm_list1', 'hmm_list2', 'transition_freq1', 'transition_freq2', 'local_div1', 'local_div2'])

	#Get .aln files
	aln_files = glob.glob(indir +'*.aln')
	
	for aln_file in aln_files:
		#Get alignment sequences and starts
		alignment, starts = parse_fasta(aln_file)
		#Match HMMs to aln
		name = aln_file.split('/')[-1]
		uids = name.split('.')[0].split('_')
		hmm_list1, null_model1, transition_freq1, local_div1 = read_hmm(indir+uids[0]+'.hhm') 
		hmm_list2, null_model2, transition_freq2, local_div2 = read_hmm(indir+uids[1]+'.hhm')
		
		#Assign sequences and starts	
		seq1 = alignment[uids[0]]
		seq2 = alignment[uids[1]]	
		s1 = int(starts[uids[0]])-1
		s2 = int(starts[uids[1]])-1

		m_hmm_list1, m_transition_freq1, m_local_div1 = match_seq_hmm(seq1, s1, hmm_list1, transition_freq1, local_div1)
		m_hmm_list2, m_transition_freq2, m_local_div2 = match_seq_hmm(seq2, s2, hmm_list2, transition_freq2, local_div2)
		pdb.set_trace()
	return None




def parse_fasta(aln_file):
	'''Parse fasta file
	'''
	alignment = {}
	starts = {}
	with open(aln_file, 'r') as file:
		for line in file:
			line = line.rstrip() #Remove newline
			if line[0] == '>':
				line = line.split('|')
				uid = line[0][1:]
				start = line[1].split(' ')[1].split('=')[1]
			else:
				alignment[uid] = line #Add sequence to dict
				starts[uid] = start
	return alignment, starts

def match_seq_hmm(seq, start, hmm_list, transition_freq, local_div):
	'''Match sequence to hmm
	'''
	m_hmm_list = []
	m_transition_freq = []
	m_local_div = []
	hmm_i = start #Which position to get from the hhm file
	for i in range(0, len(seq)):
		if seq[i] != '-':
			
			if seq[i] != hmm_list[hmm_i][0]: #If the aa don't match
				raise IOError("aa don't mtach!")

			m_hmm_list.append(hmm_list[hmm_i])
			m_transition_freq.append(transition_freq[hmm_i])
			m_local_div.append(local_div[hmm_i])
			hmm_i +=1

		else:#Gap: append zeros
			m_hmm_list.append([0]*20)
			m_transition_freq.append([0]*7)
			m_local_div.append([0]*3)
			

	return m_hmm_list, m_transition_freq, m_local_div

	
#####MAIN#####
args = parser.parse_args()
indir = args.indir[0]

#Convert alignments to HMMs and run hhblits on them
match_hh(indir)
#Parse HMMs
parse_hh(indir)
