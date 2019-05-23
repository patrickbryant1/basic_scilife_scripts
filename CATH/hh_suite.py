#! /usr/bin/env python3
# -*- coding: utf-8 -*-


import argparse
import sys
import os
import pexpect
import subprocess
import pdb

#Other script imports
from hh_reader import read_result

#FUNCTIONS
def pdb_to_fasta(uid, outdir):
	'''Convert pdb file to fasta.
	'''

	inname = uid+'.pdb'
	outname = uid+'.fa'
	#Path to pdb parser
	command = 'python /home/p/pbryant/pfs/evolution/CATH/parse_pdb_resid.py ' + inname
	outp = subprocess.check_output(command, shell = True)#Save AA sequence
	sequence = outp.split('\n')[0]
	pdb.set_trace()
	with open(outdir+outname, "w") as outfile:
		outfile.write('>'+uid+'\n')
		i = 0 #index
		while i<len(sequence):
			outfile.write(sequence[i:i+60]+'\n')
			i+=60

	return None

def run_hhblits(uid, outdir, hhblits, uniprot):
	'''A function that runs hhblits to create a HMM of an input sequence in fasta format
	'''


	inname = uid+'.fa'
	outname = uid+'.hhm'
	statement = hhblits +' -i ' +outdir+inname + ' -d ' + uniprot + ' -ohhm ' + outname
	outp = subprocess.check_output(statement, shell = True)
	
	return None

def seq_to_pdb(uids, query_aln, template_aln, start_pos, end_pos):
	'''Extracts CAs from pdb file based on sequence.
	Enables extraction of residues in alignment for
	further use.
	'''
	
	q_pdb = uids[0] + '.pdb'
	t_pdb = uids[1] + '.pdb'
	#Get query CAs
	q_command = 'python /home/p/pbryant/pfs/evolution/CATH/parse_pdb_resid.py ' + q_pdb
        q_out = subprocess.check_output(q_command, shell = True)#Save parsed pdb
	q_out = q_out.split('\n')
	q_ca = q_out[1:-1] 
	
	#Get template CAs
        t_command = 'python /home/p/pbryant/pfs/evolution/CATH/parse_pdb_resid.py ' + t_pdb
        t_out = subprocess.check_output(t_command, shell = True)#Save parsed pdb
        t_out = t_out.split('\n')
        t_ca = t_out[1:-1]   

	#Match alignment and write to file
	aligned_ca = [] #CAs that have been aligned
	for i in range(start-1, end): #Go through aligned part of sequence
		if aligned_seq[i] != '-' a§	nd ai >= (start-1) and ai < end: #If not a gap and end of alignment not reached
			aligned_ca.append(ca_pdb[ai+start]) #Append matching ca coordinates
			
	pdb.set_trace()

