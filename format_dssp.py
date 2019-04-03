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
                  default=sys.stdin, help = '''path to directory with dssp output files. The files should have names
                  on the form: uid1.phy_dssp''')

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
				'-':[0., 0., 0., 0., 0., 0., 0., 0., 0., 1.] #gap
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
			'V':165
		  }


#Functions
def parse_info(dir_path):
	'''Parse secondary structure descriptions
	and surface accessibility for dssp output.
	'''

	#chdir to directory with files
	os.chdir(dir_path)

	for file in glob.glob("*_dssp"):
		name = file.split('.')[0] #Split filename on .
		
		secondary_str = [] #save str
		surface_acc = [] #save acc

		secondary_str_hot = [] #save str in one-hot encoding




		fetch_lines = False #Don't fetch unwanted lines
		with open(file) as file:
			for line in file:
				if fetch_lines == True:
					line = line.rstrip()

					residue = line[13]
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
		

		#Format surface_acc to one_hot encodings (0-100 --> 101 possible values)
		surface_acc = np.array(surface_acc) #Convert to numpy array
		surface_acc_hot = np.eye(101)[surface_acc] #convert to one-hot encoding

		#Convert secondary_str_hot to numpy array as well
		secondary_str_hot = np.array(secondary_str_hot)

		#Write to file
		np.savetxt(name+'_acc', surface_acc_hot)
		np.savetxt(name+'_str', secondary_str_hot)

		#write_to_file(secondary_str_hot,surface_acc_hot,name)
		
	return None



def write_to_file(c1, c2, name):
	'''Write .csv file with secondary structure and surface accessibility 
	'''

	with open(name+'_dssp.csv', 'w') as file:
		for i in range(0,len(c1)):
			print(c1[i])
			pdb.set_trace()
			file.write(c1[i] + '\t' + c2[i] + '\n')



	return None


#MAIN
args = parser.parse_args()
dir_path = args.dir_path[0]

parse_info(dir_path)