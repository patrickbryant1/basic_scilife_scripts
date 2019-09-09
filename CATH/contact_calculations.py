#! /usr/bin/env python3
# -*- coding: utf-8 -*-


import argparse
import sys
import os
import glob
import numpy as np
import pandas as pd
from collections import defaultdict
from scipy.spatial import distance
import pdb

#Arguments for argparse module:
parser = argparse.ArgumentParser(description = '''Calculate contacts for both structural and sequence alignments.
						A contact is defined as 2 C-betas (C-alphas for glycine) of side-chains separated 
						by more than five residues being less than 8 Å apart.''')

parser.add_argument('indir', nargs=1, type= str,
                  default=sys.stdin, help = '''path to input directory. include / in end''')
parser.add_argument('outdir', nargs=1, type= str,
                  default=sys.stdin, help = '''path to output directory. include / in end''')
parser.add_argument('fastadir', nargs=1, type= str,
                  default=sys.stdin, help = '''path to fasta directory. include / in end''')
parser.add_argument('df_path', nargs=1, type= str,
                  default=sys.stdin, help = '''path to df.''')
parser.add_argument('hgroup', nargs=1, type= str,
                  default=sys.stdin, help = '''H-group.''')





#FUNCTIONS
def read_fasta(filepath):
	'''Read fasta sequence
	'''


	with open(filepath, 'r') as file:
		sequence = ''
		for line in file:
			line = line.rstrip() #remove \n
			if line[0] == '>':
				uid = line[1:]
			else:
				sequence += line

	return(sequence)

def match_contacts(df, indir, outdir, fastadir):
	'''Match alignments to dssp surface acc and 2ndary str
	'''


	contact_dict = {}
	fasta_dict = {}
	
	#Create new columns in df
	df['DIFFC_seqaln'] = 0
	df['DIFFC_straln'] = 0

	for index, row in df.iterrows():
		hgroup = row['H_group_x']
		uid1 = row['uid1']
		uid2 = row['uid2']
		
		if uid1 not in contact_dict.keys():
			contacts, sequence = read_cbs(indir+hgroup+'/'+uid1+'.pdb')
			contact_dict[uid1] = [contacts, sequence]
			fasta_dict[uid1] = read_fasta(fastadir+hgroup+'/'+uid1+'.fa')
		if uid2 not in contact_dict.keys():
			contacts, sequence = read_cbs(indir+hgroup+'/'+uid2+'.pdb')
			contact_dict[uid2] = [contacts, sequence]
                        fasta_dict[uid2] = read_fasta(fastadir+hgroup+'/'+uid2+'.fa')


		for suffix in ['_seqaln', '_straln']:
			
			#Match sequences to alignments
			aln1 = row['seq1'+suffix]
			aln2 = row['seq2'+suffix]
			#Save gapless alignments
			gapless_aln1 = ''
			gapless_aln2 = ''	

			for i in range(len(aln1)):
				if aln1[i] == '-' or aln2[i] == '-':
					    continue
				else:
					gapless_aln1 += aln1[i]
					gapless_aln2 += aln2[i]


			(mc1, M) = match(gapless_aln1, contact_dict[uid1], sequence_dict[uid1])
			(mc2, N) = match(gapless_aln2, contact_dict[uid2], sequence_dict[uid2])
	
			C = 0 #Keep track of the number of conserved contacts in the alignment
			for i in range(len(mc1):
				c1 = mc1[i]
				c2 = mc2[i]
				for j in c1:
					if j in c2:#If the contacts at the same position is shared.
						C+=1 
						
			diff = C/(M+N-C)
	#Write new df to outdir
	df.to_csv(outdir+hgroup+'_df.csv')
	return None

def read_cbs(pdbfile):
	'''Get the C-betas from a pdb file.
	'''
	three_to_one = {'ARG':'R', 'HIS':'H', 'LYS':'K', 'ASP':'D', 'GLU':'E', 'SER':'S', 'THR':'T', 'ASN':'N', 'GLN':'Q', 'CYS':'C', 'GLY':'G', 'PRO':'P', 'ALA':'A', 'ILE':'I', 'LEU':'L', 'MET':'M', 'PHE':'F', 'TRP':'W', 'TYR':'Y', 'VAL':'V', 'UNK': 'X'}
	sequence = ''
	pos = [] #Save positions in space
	prev_res = -1 #Keep track of potential alternative residues
	with open(pdbfile, 'r') as file:
		for line in file:
			record = parse_atm_record(line)
			if record['atm_name'] == 'CB':
				if record['res_no'] == prev_res:
					continue
				else:				
					prev_res = record['res_no'] 
					pos.append(np.array([record['x'], record['y'], record['z']]))
					sequence += three_to_one[record['res_name']] 
			if record['atm_name'] == 'CA' and record['res_name'] == 'GLY':
				prev_res = record['res_no']
				pos.append(np.array([record['x'], record['y'], record['z']]))
				sequence += three_to_one[record['res_name']]
	contacts = get_contacts(pos)
	return contacts, sequence


def parse_atm_record(line):

	record = defaultdict()
	record['name'] = line[0:6].strip()
	record['atm_no'] = int(line[6:11])
	record['atm_name'] = line[12:16].strip()
	record['res_name'] = line[17:20].strip()
	record['chain'] = line[21]
	record['res_no'] = int(line[22:26])
	record['insert'] = line[26].strip()
	record['resid'] = line[22:29]
	record['x'] = float(line[30:38])
	record['y'] = float(line[38:46])
	record['z'] = float(line[46:54])
	record['occ'] = float(line[54:60])
	record['B'] = float(line[60:66])
    
	return record

def get_contacts(pos):
	contacts = [] #Save each residue's contacts
	for i in range(len(pos)):
		contacts.append([])
		for j in range(i+5, len(pos)):
			dist = distance.euclidean(pos[i], pos[j])
			if dist < 8:
				contacts[i].append(j) 


	return contacts

def match(gapless_aln, contact_info, org_seq):
	'''Get contacts matching alignment'''


	contacts = contact_info[0]
	c_seq = contact_info[1]

	#align to original sequence
	aln1 = pairwise2.align.globalxx(org_seq, gapless_aln)
	aln2 = pairwise2.align.globalxx(org_seq, c_seq)

	seq1 = aln1[0][1] #aln to org
	seq2 = aln2[0][1] #c_seq to org

	index1 = {} #Create an index of the conversion from positions in the two sequences of the alignment
	i1=0
	for i in range(len(seq1)):
		if seq1[i] != '-':
			index1[i1]=i
			i1+=1
	
	index2 = {} #Create an index of the conversion from positions in the two sequences of the alignment
        i2=0
        for i in range(len(seq2)):
                if seq2[i] != '-':
                        index2[i2]=i
			i2+=1



	#Save matched contacts
        matched_contacts = []
	ci=0
	for i in range(len(seq1)):
		if seq1[i] != '-':
			if seq2[i] != '-': #If there are no gaps in the alignment
			else:
				matched_contacts.append([0])
		if seq2[i] != '-':
			dsspi += 1


	return matched_contacts, T

#MAIN
args = parser.parse_args()
indir = args.indir[0]
outdir = args.outdir[0]
fastadir = args.fastadir[0]
hgroup = args.hgroup[0]
df_path = args.df_path[0]
df = pd.read_csv(df_path)

df = df[df['H_group_x']==hgroup]
match_contacts(df, indir, outdir, fastadir)
