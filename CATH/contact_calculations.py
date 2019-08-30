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
parser.add_argument('df_path', nargs=1, type= str,
                  default=sys.stdin, help = '''path to df.''')





#FUNCTIONS
def match_contacts(df, indir, outdir):
	'''Match alignments to dssp surface acc and 2ndary str
	'''

	#Create empty columns in df
	df['contacts_1_seqaln'] = ''
	df['contacts_2_seqaln'] = ''
	df['contacts_1_straln'] = ''
	df['contacts_2_straln'] = ''

	contact_dict = {}
	sequence_dict = {}
	for index, row in df.iterrows():
		hgroup = row['H_group_x']
		uid1 = row['uid1']
		uid2 = row['uid2']
		
		if uid1 not in contact_dict.keys():
			contacts, sequence = read_cbs(indir+hgroup+'/'+uid1+'.pdb')
			contact_dict[uid1] = contacts
			sequence_dict[uid1] = sequence
		if uid2 not in contact_dict.keys():
			contacts, sequence = read_cbs(indir+hgroup+'/'+uid2+'.pdb')
			contact_dict[uid2] = contacts
			sequence_dict[uid2] = sequence

		for suffix in ['_seqaln', '_straln']:
			(matched_contacts) = match(row['seq1'+suffix], contact_dict[uid1], sequence_dict[uid1])
			df['contacts_1'+suffix][index] = matched_contacts			
			(matched_contacts) = match(row['seq2'+suffix], contact_dict[uid2], sequence_dict[uid2])
			df['contacts_2'+suffix][index] = matched_contacts
		pdb.set_trace()
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

def match(aln, contacts, all_seq):
	'''Get contacts matching alignment'''

	all_to_aln = {} #Save matched positions
	aln_to_all = {}
	mi = 0

	#Match all positions in aln to all_seq in order to get matching contacts
	for i in range(len(aln)):
		if mi<len(all_seq): #May be longer than sequence due to gaps
			if aln[i] == all_seq[mi]: #Matching amino acids
				
				all_to_aln[mi]=i
				aln_to_all[i]=mi
				mi += 1
		else: #If no match and the sequence is run through
			break

	#Now the positions in the full sequence is matched to the alignment
	matched_contacts = []
	for i in range(len(aln)):
		matched_contacts.append([]) #Save each residues contacts
		if i in aln_to_all.keys():#If matched - non gap
			mi = aln_to_all[i] #Get position in full seq
			for j in contacts[mi]: #Get contacts for mi
				if j in all_to_aln.keys():
					matched_contacts[i].append(all_to_aln[j])
		else:
			continue #No matching pos in all
	return matched_contacts

#MAIN
args = parser.parse_args()
indir = args.indir[0]
outdir = args.outdir[0]
df_path = args.df_path[0]
df = pd.read_csv(df_path)

match_contacts(df, indir, outdir)
