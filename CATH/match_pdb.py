#! /usr/bin/env python3
# -*- coding: utf-8 -*-


import sys
import os
#Remember to give path to singularity container
from Bio import pairwise2 
import subprocess
import pdb

def aln_to_org(org_seq, pdb_seq, aln_seq, alphas):
	'''Match fasta sequence to pdb residues
	'''
	
	failed = False
	#align to original sequence
	aln1 = pairwise2.align.globalxx(org_seq, aln_seq)
	aln2 = pairwise2.align.globalxx(org_seq, pdb_seq)
	

	if len(pdb_seq) != len(alphas):
		failed = True
	#Match aln seq to pdb
	seq1 = aln1[0][1] #aln to org
	seq2 = aln2[0][1] #pdb to org
	if len(seq1) != len(seq2):
		failed = True
	pdb_rep = ''#The matching of the pdb residues to the gapless correspondence in the alignment
	alpha_rep = []
	pdbi = 0
	for i in range(len(seq1)):
		if seq1[i] != '-': #If no gap in original aln (which has been matched to be gapless)
			if seq2[i] != '-': #If no gap in pdb aln
				pdb_rep += seq2[i]
				alpha_rep.append(alphas[pdbi])
			else:
				alpha_rep.append('-')
				pdb_rep += '-'
		if seq2[i] != '-': #If no gap in pdb aln
			pdbi += 1
			
	#org seq, aln aligned, pdb aligned			
	return pdb_rep, alpha_rep, failed

def seq_to_pdb(uids, query_aln, template_aln, q_fa, t_fa, outdir):
	'''Extracts CAs from pdb file based on sequence.
	Enables extraction of residues in alignment for
	further use.
	q_fa and q_ta are the original sequences - the full sequences for each domain.
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

	if len(q_seq) != len(q_ca):
		print(q_pdb)	
	#Get template CAs
	t_command = 'python /home/p/pbryant/pfs/evolution/CATH/parse_pdb_resid.py ' + t_pdb
	t_out = subprocess.check_output(t_command, shell = True)#Save parsed pdb
	t_out = t_out.decode()
	t_out = t_out.split('\n')
	t_seq = t_out[0]
	t_ca = t_out[1:-1]   

	if len(t_seq) != len(t_ca):
                print(t_pdb)

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
	#Match original sequence to alignment and pdb file
	q_seqmatch, q_ca_match, status_q   = aln_to_org(q_fa, q_seq, gapless_q_aln, q_ca)
	t_seqmatch, t_ca_match, status_t = aln_to_org(t_fa, t_seq, gapless_t_aln, t_ca)	
	
	if status_q == False and status_t == False:
		#Match alignment and write to file
		q_file = open(outdir+uids[0]+'_to_'+uids[1]+'_aln.pdb', 'w')
		t_file = open(outdir+uids[1]+'_to_'+uids[0]+'_aln.pdb', 'w')
	
		#Keep track of true index
		ti = 0

		for i in range(0, len(q_seqmatch)): #Go through aligned part of sequence and only select residues when both sequencenes do not have a gap in extracted pdb alignment
			write_to_file = False #Keep track of if or not to write to at each position
			if q_seqmatch[i] != '-' and t_seqmatch[i] != '-': #No gap in either query or template
				write_to_file = True
				ti+=1
			if write_to_file == True:
				replace_str = ' '+str(ti)+' '*(8-len(str(ti)))
				q_file.write(q_ca_match[i][0:22]+replace_str+q_ca_match[i][31:]+'\n') #Write matching ca coordinates - the ca and seq should be matched
				t_file.write(t_ca_match[i][0:22]+replace_str+t_ca_match[i][31:]+'\n')

		q_file.close()
		t_file.close()
	return None
