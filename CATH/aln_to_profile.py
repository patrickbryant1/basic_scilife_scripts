#! /usr/bin/env python3
# -*- coding: utf-8 -*-


import argparse
import sys
import glob

#Custom import
from parse_hmm import read_hmm



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


def parsr_hh(indir):
	'''Parse HMMs
	'''

	#Get aligned hmm files blitsed against uiprot
        files = glob.glob(input_dir +'*re.hhm')

	for file in files:
		read_hmm
		
#####MAIN#####
args = parser.parse_args()

indir = args.indir[0]
hhmake = args.hhmake[0]
hhblits = args.hhblits[0]
uniprot = args.uniprot[0]

#Convert alignments to HMMs
run_hh(indir, hhmake, hhblits, uniprot)
