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
								each residue.''')
 
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
		uids = name[0].split('_') #Separate uids
		uid_1 = uids[0] #Get uid_1
		uid_2 = uids[1] #Get uid_2
		
		with open(file) as file:
			for line in file:
				line = line.rstrip()









#MAIN
args = parser.parse_args()
dir_path = args.dir_path[0]
