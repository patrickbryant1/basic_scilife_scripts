#! /usr/bin/env python2
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

	inname = outdir+uid+'.pdb'
	outname = inname+'.fa'
	#Path to pdb parser
	command = 'python /home/p/pbryant/pfs/evolution/CATH/parse_pdb_resid.py ' + inname
	outp = subprocess.check_output(command, shell = True)#Save AA sequence
	outp = outp.decode()
	sequence = outp.split('\n')[0]
	#pdb.set_trace()
	with open(outname, "w") as outfile:
		outfile.write('>'+uid+'\n')
		i = 0 #index
		while i<len(sequence):
			outfile.write(sequence[i:i+60]+'\n')
			i+=60

	return None

def run_hhblits(uid, indir, hhblits, uniprot):
	'''A function that runs hhblits to create a HMM of an input sequence in fasta format
	'''


	inname = indir+uid+'.fa'
	outname = uid+'.hhm'
	statement = hhblits +' -i ' + inname + ' -d ' + uniprot + ' -ohhm ' + outname
	outp = subprocess.check_output(statement, shell = True)
	
	return None
def match_aln_pdb(pdb_seq, aln_seq, start, end):
	'''Match fasta sequence to pdb residues
	'''
	
	pdb_rep = '' #Save pdb representation
	#Go through alignment and create pdb representation of the same alignment
	pdb_i = start-1
	for i in range(0, len(aln_seq)):
		if aln_seq[i] == pdb_seq[pdb_i]: #If match - add
			pdb_rep += aln_seq[i]
			pdb_i+=1
		else:
			pdb_rep += '-' #If no match - add gap
 
		print(aln_seq[i], pdb_seq[pdb_i], pdb_rep)
	pdb.set_trace()
	return pdb_rep


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
	q_out = q_out.decode() #Returns byte
	q_out = q_out.split('\n')
	q_seq = q_out[0]	
	q_ca = q_out[1:-1] 
	
	#Get template CAs
	t_command = 'python /home/p/pbryant/pfs/evolution/CATH/parse_pdb_resid.py ' + t_pdb
	t_out = subprocess.check_output(t_command, shell = True)#Save parsed pdb
	t_out = t_out.decode()
	t_out = t_out.split('\n')
	t_seq = t_out[0]
	t_ca = t_out[1:-1]   
	
	pdb.set_trace()
	#Create representation of alignment due to pdb file
	q_match = match_aln_pdb(q_seq, query_aln, start_pos[0], end_pos[0])

	#Match alignment and write to file
	q_file = open(uids[0]+'to'+uids[1]+'_aln.pdb', 'w')
	t_file = open(uids[1]+'to'+uids[0]+'_aln.pdb', 'w')
	q_start = start_pos[0]
	t_start = start_pos[1]
	
 
	qi = q_start-2#residue to get
	ti = t_start-2
	pdb.set_trace()

	
	for i in range(0, len(query_aln)): #Go through aligned part of sequence and only select residues when both sequencenes do not have a gap
		write_to_file = False #Keep track of if or not to write to at each position
		if query_aln[i] != '-' and template_aln[i] != '-': #No gap in either query or template
			qi+=1 #Increase indexes, since no gaps
			ti+=1
			write_to_file = True
		if template_aln[i] == '-': #If a gap in template
			qi+=1 #Increase query index
		if query_aln[i] == '-': #If a gap in query
                        ti+=1 #Increase template index
		if write_to_file == True:
			if q_seq[qi] == query_aln[i]:#Check AAs match
				q_file.write(q_ca[qi]+'\n') #Write matching ca coordinates
			else:
				raise ValueError('query pos: ' + str(qi))

			if t_seq[ti] == template_aln[i]:
				t_file.write(t_ca[ti]+'\n')
			else:
				raise ValueError('template pos: ' + str(ti))

	q_file.close()
	t_file.close()
	return None

def make_phylip(uids, query_aln, template_aln):
        '''Print phylip format for tree-puzzle calculations
        '''
        #Create text in phylip format
        text = (' 4  ' + str(len(query_aln)) + '\n'
                        + uids[0] + '00|' + query_aln + '\n'
                        + 'copy11111' + '|' + query_aln + '\n'
                        + uids[1] + '00|' + template_aln + '\n'
                        + 'copy22222' + '|' + template_aln + '\n')


        #Define file name
        file_name = uids[0] + '_' + uids[1] + '.phy'
        #Open file and write text to it
        with open(file_name, "w") as file:
                file.write(text)

        return None

