#! /usr/bin/env python3
# -*- coding: utf-8 -*-


import argparse
import sys
import os
import subprocess
import pdb

#Other script imports
from hh_reader import read_result


#FUNCTIONS
def pdb_to_fasta(uid, outdir):
	'''Convert pdb file to fasta.
	'''
	letters = {'ALA':'A','ARG':'R','ASN':'N','ASP':'D','CYS':'C','GLU':'E','GLN':'Q','GLY':'G','HIS':'H',
           	'ILE':'I','LEU':'L','LYS':'K','MET':'M','PHE':'F','PRO':'P','SER':'S','THR':'T','TRP':'W',
		'TYR':'Y','VAL':'V'}

	inname = uid+'.pdb'
	outname = uid+'.fa'

	sequence = '' #Save AA sequence
	with open(inname, "r") as infile:
		for line in infile:
			line = line.split()
			if line[0] == 'ATOM' and line[2] == 'CA':
				sequence+=line[3]

	with open(outdir+outname, "w") as outfile:
		outfile.write('>'+uid+'\n')
		i = 0 #index
		while i<len(sequence):
			outfile.write(sequence[i:i+60]+'\n')
			i+=60

	return None
