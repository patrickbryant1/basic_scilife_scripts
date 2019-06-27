#! /usr/bin/env python3
# -*- coding: utf-8 -*-


import argparse
import sys
import glob
import pandas as pd
#Custom import
from parse_hmm import read_hmm
#Debugging
import pdb


#Arguments for argparse module:
parser = argparse.ArgumentParser(description = '''A program that converts a fasta alignment to a HMM and runs hhblits on it.
						The output from hhblits is parsed and written to a .csv file.''')

parser.add_argument('indir', nargs=1, type= str,
default=sys.stdin, help = 'path to directory with .aln files.')

parser.add_argument('hhmake', nargs=1, type= str,
default=sys.stdin, help = 'path to hhmake.')

parser.add_argument('hhblits', nargs=1, type= str,
default=sys.stdin, help = 'path to hhblits.')

parser.add_argument('uniprot', nargs=1, type= str,
default=sys.stdin, help = 'path to uniprot database.')

#FUNCTIONS
def run_hh(indir, hhmake, hhblits, uniprot):
	'''Convert alignments to HMMs, then run
	hhblits on HMM and parse output to csv'''

	#Get .aln files
	aln_files = glob.glob(input_dir +'*.aln')

	for file in aln_files:
		#Make HMM of aln
		name = file.split('/')[-1]
		call_hhmake = hhmake +' -i ' + file + ' -M 50 ' + ' -o ' + name+'.hhm'
		outp = subprocess.check_output(call_hhmake, shell = True)
		
		#Run hhblits on new HMM representation of aln
		call_hhblits = hhblits +' -i ' + name + '.hhm -d ' + uniprot + ' -ohhm ' + name + '.re.hhm'
		outp = subprocess.check_output(call_hhblits, shell = True)
		pdb.set_trace()
	return None

def parse_hh(indir):
	'''Parse HMMs
	'''
	#Pandas df
	df = pd.DataFrame(columns = ['uid1','uid2','hmm_list', 'null_model', 'transition_freq', 'local_div'])
	#Get aligned hmm files blitsed against uiprot
        files = glob.glob(input_dir +'*re.hhm')
	

	for file in files:
		hmm_list, null_model, transition_freq, local_div = read_hmm(file)
		new_df = pd.DataFrame([hmm_list, null_model, transition_freq[1:], local_div[1:]], columns = ['uid1','uid2','hmm_list', 'null_model', 'transition_freq', 'local_div'])
		df.append([hmm_list, null_model, transition_freq[1:], local_div[1:]], ignore_index = True)


	df.to_csv('hmm_df.csv',index=False)

	return None
		
#####MAIN#####
args = parser.parse_args()

indir = args.indir[0]
hhmake = args.hhmake[0]
hhblits = args.hhblits[0]
uniprot = args.uniprot[0]

#Convert alignments to HMMs and run hhblits on them
run_hh(indir, hhmake, hhblits, uniprot)
#Parse HMMs
parse_hh(indir)
