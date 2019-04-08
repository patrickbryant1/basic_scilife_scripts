#! /usr/bin/env python3
# -*- coding: utf-8 -*-


import argparse
import sys
import os
import glob
import numpy as np


import pdb


#Arguments for argparse module:
parser = argparse.ArgumentParser(description = '''A program that parses output from dssp and returns 
								secondary structure descriptions accordingly:
								    G = 3-turn helix (310 helix). Min length 3 residues.
    								H = 4-turn helix (α helix). Minimum length 4 residues.
    								I = 5-turn helix (π helix). Minimum length 5 residues.
    								T = hydrogen bonded turn (3, 4 or 5 turn)
   									E = extended strand in parallel and/or anti-parallel β-sheet conformation. Min length 2 residues.
    								B = residue in isolated β-bridge (single pair β-sheet hydrogen bond formation)
    								S = bend (the only non-hydrogen-bond based assignment).
    								C = coil (residues which are not in any of the above conformations). 
    							
    							in a one-hot encoding where each of the eight secondary structure elements as well as missing secondary
								structure predictions and gaps in alignments are represented. E.g. G = [1000000000].

								Surface accessibility for each residue is also reported in one-hot encodings in the range of 0-100 percent
								after being normalized due to the maximum accessible surface areas for each amino acid according to empirical 
								measurements in: Tien, Matthew Z et al. “Maximum allowed solvent accessibilites of residues in proteins.” 
								PloS one vol. 8,11 e80635. 21 Nov. 2013, doi:10.1371/journal.pone.0080635. 
								''')
 
parser.add_argument('dir_path', nargs=1, type= str,
                  default=sys.stdin, help = '''path to directory with dssp output and TMalign per residue alignment files. 
                  The dssp files should have names on the form: uid1.phy_dssp
                  The alignment files should have names on the form: uid1_uid2.aln ''')

parser.add_argument('out_path', nargs=1, type= str,
                  default=sys.stdin, help = '''path to output directory. include / in end''')



####Normalization and Encoding####

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

#Functions
def encode_dssp(dir_path):
	'''Parse secondary structure descriptions
	and surface accessibility for dssp output.
	'''

	#chdir to directory with files
	os.chdir(dir_path)

	dssp_info = {} #save dssp output for encodings

	for file in glob.glob("*_dssp"):
		uid = file.split('.')[0] #Split filename on .
		
		secondary_str = [] #save str
		surface_acc = [] #save acc




		fetch_lines = False #Don't fetch unwanted lines
		with open(file) as file:
			for line in file:
				if fetch_lines == True:
					
					line = line.rstrip()
					str_i = line[16]
                    residue = line[13]
					acc_i = line[35:38].strip()
					
					if residue != '!' and residue != '!*':
					                                           
						#Normalize acc_i by the max acc surface area for the specific amino acid
						#Round to whole percent
						acc_i_norm = round((float(acc_i)/max_acc[residue])*100, )
						acc_i_norm = min(acc_i_norm, 100) #Should not be over 100 percent
					
					else: #If a chain break
						acc_i_norm = 0
					

					#Add values
					secondary_str.append(str_i)
					surface_acc.append(acc_i_norm)

				if '#' in line:
					fetch_lines = True
					#now the subsequent lines will be fetched
		

		#Save to dict
		dssp_info[uid] = [secondary_str, surface_acc]
		
		
	return dssp_info


def match_encoding(seq1, dssp1, seq2, dssp2):
	'''Make encoding of amino acids in sequence and add corresponding encoding for 
	dssp metrics and write to file.
	'''

	encoded_aln = [] #Save the complete encoding
	#assign structrural encodings and normalized acc
	str1, acc1 = dssp1[0], dssp1[1]
	str2, acc2 = dssp2[0], dssp2[1]

	pos1 = 0 #keep track of dssp pos, will not match to sequence due to gaps
	pos2 = 0

	gap_acc = 0 #The gap surface accessibility should be 0
	enc = [] #save all encodings for each residue

	print(len(seq1))

	for i in range(0,len(seq1)):
		aa1 = seq1[i] #Get residues
		aa2 = seq2[i]
		print(pos1,len(str1), pos2, len(str2))
		if pos2 >= len(str2):
			pdb.set_trace()
		if aa1 == '-': #if a gap in 1
			aa1_str = '-' #set encoding to gap
			aa1_acc = gap_acc #get gap acc

			enc = [aa1, aa1_str, gap_acc, aa2, str2[pos2], acc2[pos2]] #Encode as pair (order is important! Otherwise the order of the subsequent residues will not be preserved)
			encoded_aln.append(enc) #Append encoding to full alignment encoding
			pos2 +=1 #Increase pos 2 (pos1 is gap)

		elif aa2 == '-': #if a gap in 2
			aa2_str = '-' #set encoding to gap
			aa2_acc = gap_acc #get gap acc

			enc = [aa1, str1[pos1], acc1[pos1], aa2, aa2_str, aa2_acc] #Encode as pair (order is important! Otherwise the order of the subsequent residues will not be preserved)
			encoded_aln.append(enc) #Append encoding to full alignment encoding
			pos1 +=1 #Increase pos 2 (pos1 is gap)
		
		else: #if something else than a gap in both
			enc = [aa1, str1[pos1], acc1[pos1], aa2, str2[pos2], acc2[pos2]] #Encode as pair (order is important! Otherwise the order of the subsequent residues will not be preserved)
			encoded_aln.append(enc) #Append encoding to full alignment encoding
			pos1 +=1
			pos2 +=1


	return encoded_aln

def encode_aln(dir_path, dssp_info, out_path):

	for file in glob.glob("*.aln"):
		names = file.split('.')[0] #Split filename on .
		names = names.split('_')
		uid1 = names[0]
		uid2 = names[1]

		with open(file) as file:
			n = 0 #keep track of line number
			for line in file:
				if n == 0:
					seq1 = line.rstrip()
					n+=1
				if n == 1:
					seq2 = line.rstrip()


		if len(seq1) == len(seq2): #If the seqs are of equal lengths
			#Get dssp encodings for uids
			dssp1 = dssp_info[uid1]
			dssp2 = dssp_info[uid2]
		else:
			raise ValueError('Alignments are of different lengths for: ' + uid1, uid2)


		

		print(uid1, uid2) #Print uid pairs
		encoded_aln = match_encoding(seq1, dssp1, seq2, dssp2)

		name = out_path+uid1+'_'+uid2+'.enc'
		write_encoding(encoded_aln, name)

	return None



def write_encoding(encoded_aln, name):
	'''Write the encoding to a file
	'''

	with open(name, 'w') as file:
		for enc_pair in encoded_aln:
			for i in enc_pair:
				file.write(str(i)+',')
			file.write('\n')


	return None
#MAIN
args = parser.parse_args()
dir_path = args.dir_path[0]
out_path = args.out_path[0]

#Encode dssp output
dssp_info = encode_dssp(dir_path)

#Encode alignment and add dssp encoding and write to file
encode_aln(dir_path, dssp_info, out_path)

