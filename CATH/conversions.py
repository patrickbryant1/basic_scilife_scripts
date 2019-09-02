#! /usr/bin/env python2
# -*- coding: utf-8 -*-


import argparse
import sys
import os
import pexpect
import subprocess
#Remember to give path to singularity container
from Bio import pairwise2 
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

def match_aln_pdb(pdb_seq, alphas, aln_seq):
	'''Match fasta sequence to pdb residues
	'''
	
	pdb_rep = '' #Save pdb representation (the ones that match)
	alpha_rep = [] #Save represented CAs (the ones that match)
	#Go through alignment and create pdb representation of the same alignment

	#align
	alignments = pairwise2.align.globalxx(aln_seq, pdb_seq)
	seq1 = alignments[0][0] #Will be the alignment of aln_seq to pdb_seq
	seq2 = alignments[0][1] #Will be the alignment of pdb_seq to aln_seq

	#pdb index
	pdb_i = 0
	for i in range(0, len(seq1)):
		if seq1[i] != '-' and seq2[i] != '-':
			pdb_rep += seq2[i]
			alpha_rep.append(alphas[pdb_i])
		if seq2[i] != '-': #If no gap - increase pdb index
                        pdb_i += 1
			
	return pdb_rep, alpha_rep

def seq_to_pdb(uids, query_aln, template_aln, start_pos, end_pos, outdir):
	'''Extracts CAs from pdb file based on sequence.
	Enables extraction of residues in alignment for
	further use.
	'''
	
	q_pdb = uids[0] + '.pdb'
	t_pdb = uids[1] + '.pdb'
	
	#Get the .pdb if they do not exist
	for uid in uids:
		if not os.path.isfile(uid+'.pdb'):
			subprocess.call(["wget",'www.cathdb.info/version/v4_2_0/api/rest/id/'+uid+'.pdb'])
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

	#Save gapless alignments
	gapless_q_aln = ''
	gapless_t_aln = ''
	for i in range(0, len(query_aln)): #Go through aligned part of sequence and only select residues when both sequencenes do not have a gap in extracted alignment
		if query_aln[i] == '-' or template_aln[i] == '-' :
			continue
		else:
			gapless_q_aln += query_aln[i]
			gapless_t_aln += template_aln[i]

	
	#Create representation of gapless alignments due to pdb file
	q_seq_match, q_ca_match = match_aln_pdb(q_seq, q_ca, gapless_q_aln)
	t_seq_match, t_ca_match = match_aln_pdb(t_seq, t_ca, gapless_t_aln)	
	#Match alignment and write to file
	q_file = open(outdir+uids[0]+'_to_'+uids[1]+'_aln.pdb', 'w')
	t_file = open(outdir+uids[1]+'_to_'+uids[0]+'_aln.pdb', 'w')
	
	pdb.set_trace()	
	for i in range(0, len(q_seq_match)): #Go through aligned part of sequence and only select residues when both sequencenes do not have a gap in extracted pdb alignment
		write_to_file = False #Keep track of if or not to write to at each position
		if q_seq_match[i] != '-' and t_seq_match[i] != '-': #No gap in either query or template
			write_to_file = True
		if write_to_file == True:
			replace_str = ' '+str(i)+' '*(8-len(str(i)))
			q_file.write(q_ca_match[i][0:22]+replace_str+q_ca_match[i][31:]+'\n') #Write matching ca coordinates
			t_file.write(t_ca_match[i][0:22]+replace_str+t_ca_match[i][31:]+'\n')

	q_file.close()
	t_file.close()
	return None

def make_phylip(uids, query_aln, template_aln, outdir):
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
        with open(outdir+file_name, "w") as file:
                file.write(text)

        return None

