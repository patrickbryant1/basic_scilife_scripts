#! /usr/bin/env python3
# -*- coding: utf-8 -*-


import argparse
import sys
import os
import glob
import numpy as np
import pandas as pd
import pdb


#Arguments for argparse module:
parser = argparse.ArgumentParser(description = '''Parse dssp and match to both structural and sequence alignments''')

parser.add_argument('indir', nargs=1, type= str,
                  default=sys.stdin, help = '''path to input directory. include / in end''')
parser.add_argument('outdir', nargs=1, type= str,
                  default=sys.stdin, help = '''path to output directory. include / in end''')
parser.add_argument('df_path', nargs=1, type= str,
                  default=sys.stdin, help = '''path to df.''')

#Max acc surface areas for each amino acid according to empirical measurements in:
#Tien, Matthew Z et al. “Maximum allowed solvent accessibilites of residues in proteins.” 
#PloS one vol. 8,11 e80635. 21 Nov. 2013, doi:10.1371/journal.pone.0080635
max_acc = { 'A':121,
			'R':265,
			'N':187,
			'D':187,
			'C':148,
			'E':214,
			'Q':214,
			'G':97,
			'H':216,
			'I':195,
			'L':191,
			'K':230,
			'M':203,
			'F':228,
			'P':154,
			'S':143,
			'T':163,
			'W':264,
			'Y':255,
			'V':165,
			'X':192 #Average of all other maximum surface accessibilites
		  }



#FUNCTIONS
def match_dssp_to_aln(df, indir, outdir):
	'''Match alignments to dssp surface acc and 2ndary str
	'''

	dssp_dict = {}

	#Create new columns in df
	df['2ndarystr_1_seqaln'] = ''
	df['acc_1_seqaln'] = ''
	df['2ndarystr_2_seqaln'] = ''
	df['acc_2_seqaln'] = ''
	df['2ndarystr_1_straln'] = ''
	df['acc_1_straln'] = ''
	df['2ndarystr_2_straln'] = ''
	df['acc_2_straln'] = ''

	for index, row in df.iterrows():
		hgroup = row['H_group_x']
		uid1 = row['uid1']
		uid2 = row['uid2']

		if uid1 not in dssp_dict.keys():
			sequence1, secondary_str1, surface_acc1 = parse_dssp(indir+hgroup+'/dssp/'+uid1+'.pdb.dssp')
			dssp_dict[uid1] = [sequence1, secondary_str1, surface_acc1]
		if uid2 not in dssp_dict.keys():
			sequence2, secondary_str2, surface_acc2 = parse_dssp(indir+hgroup+'/dssp/'+uid2+'.pdb.dssp')
			dssp_dict[uid2] = [sequence2, secondary_str2, surface_acc2]

		for suffix in ['_seqaln', '_straln']:
			if suffix == '_seqaln':
				start1 = row['s1']
				start2 = row['s2']
			else:
				start1 = 0
				start2 = 0
			(matched_secondary_str, matched_surface_acc) = match(row['seq1'+suffix], dssp_dict[uid1], start1)
			df['2ndarystr_1'+suffix][index] = matched_secondary_str
			df['acc_1'+suffix][index] = matched_surface_acc

			(matched_secondary_str, matched_surface_acc) = match(row['seq2'+suffix], dssp_dict[uid2], start2)
			df['2ndarystr_2'+suffix][index] = matched_secondary_str
			df['acc_2'+suffix][index] = matched_surface_acc

	#Write new df to outdir
	df.to_csv(outdir+'complete_df.csv')
	return None

def parse_dssp(filepath):
	'''Parse dssp output and return sequence, secondary structure and surface accessibilities.
	'''

	secondary_str = [] #save str
	surface_acc = [] #save acc
	sequence = '' #save residues
	fetch_lines = False #Don't fetch unwanted lines
	try:
		with open(filepath, 'r') as file:
			for line in file:
				if fetch_lines == True:
					line = line.rstrip()
					residue = line[13]
					str_i = line[16]
					acc_i = line[35:38].strip()
				
					if residue == '!' or residue == '!*':
						continue
					else:	                                           
					#Normalize acc_i by the max acc surface area for the specific amino acid
					#Round to whole percent
						acc_i_norm = round((float(acc_i)/max_acc[residue])*100, )
						acc_i_norm = min(acc_i_norm, 100) #Should not be over 100 percent
					
						#Save
						secondary_str.append(str_i)
						surface_acc.append(acc_i_norm)
						sequence = sequence + residue
				if '#' in line:
					fetch_lines = True
					#now the subsequent lines will be fetched
	except:
		if os.path.isfile(filepath):
			raise ValueError('File', filepath,'exists but could not be read')

	return(sequence, secondary_str, surface_acc)						

def match(aln, info, start):
	'''Extract secondary structure annotations and surface acc matching the aligned sequences'''

	all_seq = info[0]
	all_str = info[1]
	all_acc = info[2]

	matched_str = ''
	matched_acc = ''

	mi = start-1
	print(aln[0], all_seq[mi])
	for i in range(len(aln)):
		if mi<len(all_seq): #May be longer than sequence due to gaps
			if aln[i] == all_seq[mi]:
				matched_str += all_str[mi]+','	
				matched_acc += str(all_acc[mi])+','
				mi += 1
			else:
				matched_str += '-'+','
				matched_acc += '-'+','
		else: #If no match and the sequence is run through - add gap
			matched_str += '-'+','
			matched_acc += '-'+','
	return (matched_str, matched_acc)

#MAIN
args = parser.parse_args()
indir = args.indir[0]
outdir = args.outdir[0]
df_path = args.df_path[0]
df = pd.read_csv(df_path)

match_dssp_to_aln(df, indir, outdir)

