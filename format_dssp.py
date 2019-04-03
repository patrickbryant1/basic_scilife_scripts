#! /usr/bin/env python3
# -*- coding: utf-8 -*-


import argparse
import sys
import os
import glob
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
    							As well as surface accessibility for
								each residue. The output is in .csv format''')
 
parser.add_argument('dir_path', nargs=1, type= str,
                  default=sys.stdin, help = '''path to directory with dssp output files. The files should have names
                  on the form: uid1.phy_dssp''')


#Functions
def parse_info(dir_path):
	'''Parse secondary structure descriptions
	and surface accessibility for dssp output.
	'''

	#chdir to directory with files
	os.chdir(dir_path)

	for file in glob.glob("*_dssp"):
		name = file.split('.')[0] #Split filename on .
		
		c1 = [] #column 1 in output
		c2 = [] #column 2 in output

		fetch_lines = False #Don't fetch unwanted lines
		with open(file) as file:
			for line in file:
				if fetch_lines == True:
					line = line.rstrip()

					secondary_str = line[16]
					surface_acc = line[35:38].strip()

					c1.append(secondary_str)
					c2.append(surface_acc)
					

				if '#' in line:
					fetch_lines = True
					#now the subsequent lines will be fetched
		write_to_file(c1,c2,name)


	return None

def write_to_file(c1, c2, name):
	'''Write .csv file with secondary structure and surface accessibility 
	'''

	with open(name+'_dssp.csv', 'w') as file:
		for i in range(0,len(c1)):
			file.write(c1[i] + ',' + c2[i] + '\n')



	return None


#MAIN
args = parser.parse_args()
dir_path = args.dir_path[0]

parse_info(dir_path)