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

parser.add_argument('df_path', nargs=1, type= str,
                  default=sys.stdin, help = '''path to df.''')

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
def match_dssp_to_aln(df, indir):
	'''Match alignments to dssp surface acc and 2ndary str
	'''

	for index, row in df.iterrows():
		hgroup = row['H_group._x']
		uid1 = row['uid1']
		uid2 = row['uid2']

		sequence1, secondary_str1, surface_acc1 = parse_dssp(indir+hgroup+'/dssp/'+uid1+'.dssp')
		sequence2, secondary_str2, surface_acc2 = parse_dssp(indir+hgroup+'/dssp/'+uid2+'.dssp')
#MAIN
args = parser.parse_args()
indir = args.indir[0]
df_path = args.df_path[0]
df = pd.read_csv(df_path)



