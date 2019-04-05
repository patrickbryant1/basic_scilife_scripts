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
#one-hot encoding
str_encoding = {'G':[1., 0., 0., 0., 0., 0., 0., 0., 0., 0.],
				'H':[0., 1., 0., 0., 0., 0., 0., 0., 0., 0.],
				'I':[0., 0., 1., 0., 0., 0., 0., 0., 0., 0.],
				'T':[0., 0., 0., 1., 0., 0., 0., 0., 0., 0.],
				'E':[0., 0., 0., 0., 1., 0., 0., 0., 0., 0.],
				'B':[0., 0., 0., 0., 0., 1., 0., 0., 0., 0.],
				'S':[0., 0., 0., 0., 0., 0., 1., 0., 0., 0.],
				'C':[0., 0., 0., 0., 0., 0., 0., 1., 0., 0.],
				' ':[0., 0., 0., 0., 0., 0., 0., 0., 1., 0.],
				'-':[0., 0., 0., 0., 0., 0., 0., 0., 0., 1.] #gap, used in later stage
				}

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

#Encoding for amino acid sequence from extracted per residue alignment from TMalign structural alignment
seq_encoding = {'A':[1., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0.],
				'R':[0., 1., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0.],
				'N':[0., 0., 1., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0.],
				'D':[0., 0., 0., 1., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0.],
				'C':[0., 0., 0., 0., 1., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0.],
				'E':[0., 0., 0., 0., 0., 1., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0.],
				'Q':[0., 0., 0., 0., 0., 0., 1., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0.],
				'G':[0., 0., 0., 0., 0., 0., 0., 1., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0.],
				'H':[0., 0., 0., 0., 0., 0., 0., 0., 1., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0.],
				'I':[0., 0., 0., 0., 0., 0., 0., 0., 0., 1., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0.],
				'L':[0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 1., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0.],
				'K':[0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 1., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0.],
				'M':[0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 1., 0., 0., 0., 0., 0., 0., 0., 0., 0.],
				'F':[0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 1., 0., 0., 0., 0., 0., 0., 0., 0.],
				'P':[0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 1., 0., 0., 0., 0., 0., 0., 0.],
				'S':[0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 1., 0., 0., 0., 0., 0., 0.],
				'T':[0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 1., 0., 0., 0., 0., 0.],
				'W':[0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 1., 0., 0., 0., 0.],
				'Y':[0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 1., 0., 0., 0.],
				'V':[0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 1., 0., 0.],
				'X':[0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 1., 0.],
				'-':[0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 1.] 
				}
#Functions
def encode_dssp(dir_path):
	'''Parse secondary structure descriptions
	and surface accessibility for dssp output.
	'''

	#chdir to directory with files
	os.chdir(dir_path)

	dssp_hot = {} #save dssp output in one-hot encodings

	for file in glob.glob("*_dssp"):
		uid = file.split('.')[0] #Split filename on .
		
		secondary_str = [] #save str
		surface_acc = [] #save acc

		secondary_str_hot = [] #save str in one-hot encoding




		fetch_lines = False #Don't fetch unwanted lines
		with open(file) as file:
			for line in file:
				if fetch_lines == True:
					line = line.rstrip()

					residue = line[13]
					if residue != '!' and residue != '!*':
						str_i = line[16]
						acc_i = line[35:38].strip()

						#Normalize acc_i by the max acc surface area for the specific amino acid
						#Round to whole percent
						acc_i_norm = round((float(acc_i)/max_acc[residue])*100, )
						acc_i_norm = min(acc_i_norm, 100) #Should not be over 100 percent
						#Add original values
						secondary_str.append(str_i)
						surface_acc.append(acc_i_norm)

						#Add one-hot encodings
						secondary_str_hot.append(str_encoding[str_i])
				
					

				if '#' in line:
					fetch_lines = True
					#now the subsequent lines will be fetched
		

		#Format surface_acc to one_hot encodings (0-100, gap, unknown residues. The
		#gap and unknown residues will be assigned 0 --> 101 possible values)
		surface_acc = np.array(surface_acc) #Convert to numpy array
		surface_acc_hot = np.eye(101)[surface_acc] #convert to one-hot encoding

		#Convert secondary_str_hot to numpy array as well
		secondary_str_hot = np.array(secondary_str_hot)


		#Save to dict
		dssp_hot[uid] = [surface_acc_hot, secondary_str_hot]
		
		
	return dssp_hot


def encode_aln(dir_path, dssp_hot, out_path):

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

		#Get one-hot dssp encodings for uids
		dssp1 = dssp_hot[uid1]
		dssp2 = dssp_hot[uid2]

		print(uid1, uid2)
		enc1 = encode_out(seq1, dssp1)
		enc2 = encode_out(seq2, dssp2)

		if len(enc1) == len(enc2): #If the encodings are of equal lengths
			aln_matrix = []

			for i in range(0, len(enc1)):
				aln_matrix.append(np.concatenate((enc1[i], enc2[i])))


		else:
			raise ValueError('Encodings are of different lengths for: ' + uid1, uid2)


		aln_array = np.array(aln_matrix)

		#Save matrix to disk
		
		np.savetxt(out_path+uid1+'_'+uid2+'.hot', aln_array)#, fmt='%d')
		#pdb.set_trace()
	return None

def encode_out(sequence, dssp):
	'''Make one-hot encoding of amino acids in sequence
	and add corresponding encoding for dssp metrics and
	write to file
	'''

	#assign acc and structrural encodings
	surface_acc_hot, secondary_str_hot = dssp[0], dssp[1]
	dssp_a = 0 #keep track of dssp pos, will not match to sequence due to gaps

	gap_acc = np.eye(101)[0] #The gap surface accessibility should be 0
	all_hot = [] #save all one-hot encodings for each residue


	for a in sequence:
		a_hot = seq_encoding[a] #Encode residue

		if a == '-': #if a gap
			a_str = str_encoding[a] #get str encoding for gap
			a_acc = gap_acc #get gap acc
			cat_enc = np.concatenate((a_hot, a_str, a_acc)) #cat
			all_hot.append(cat_enc) #append to representation

		else: #if something else than a gap
			a_str = secondary_str_hot[dssp_a] #get str encoding 
			a_acc = surface_acc_hot[dssp_a] #get acc encoding
			cat_enc = np.concatenate((a_hot, a_str, a_acc)) #cat
			all_hot.append(cat_enc) #append to representation
			

			dssp_a +=1

	return all_hot

def write_encoding(aln_matrix, name):
	'''Write the one-hot encoding to a file
	'''

	with open(name, 'w') as file:
		for i in aln_matrix:
			file.write(i)

#MAIN
args = parser.parse_args()
dir_path = args.dir_path[0]
out_path = args.out_path[0]

#Encode dssp output
dssp_hot = encode_dssp(dir_path)

#Encode alignment and add dssp encoding and write to file
encode_aln(dir_path, dssp_hot, out_path)

