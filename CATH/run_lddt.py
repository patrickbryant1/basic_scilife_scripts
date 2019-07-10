#! /usr/bin/env python3
# -*- coding: utf-8 -*-


import argparse
import sys
import os
import subprocess
import glob
import pdb

#Arguments for argparse module:
parser = argparse.ArgumentParser(description = '''A program that renumbers the entrties in pdb files
						and then runs lddt on the pdb files.''')

parser.add_argument('indir', nargs=1, type= str,
default=sys.stdin, help = 'path to directory with pdb files.')


def run_lddt(indir):
	'''Renumber pdb files to start with atom 0, residue 1.
	'''

	pdb_files = glob.glob(indir +'*_aln.pdb')
	
	#Run LDDT
	while pdb_files:
		pdb_name1 = pdb_files[0].split('/')[-1]
		uid1 = pdb_name1.split('_')[0]
		uid2 = pdb_name1.split('_')[2]
		
		for i in range(1, len(pdb_files)):
			if uid1 in pdb_files[i] and uid2 in pdb_files[i]:
				break
		pdb_name2 = pdb_files[i].split('/')[-1]
		#Fix pdb files
		move_res_number(pdb_name1)
		move_res_number(pdb_name2)
		#Run lddt on fixed pdb files
		command = 'singularity run --app lDDT /home/p/pbryant/pfs/singularity/ost.img -c -x -t ' + 'rf_'+pdb_name1 + ' rf_'+pdb_name2
		out = subprocess.check_output(command, shell = True)#Save parsed pdb
		
		#Write to file
		write_lddt(indir, out, uid1, uid2)

		#Pop used pdb files
		pdb_files.pop(i)
		pdb_files.pop(0)
		

	return None

def move_res_number(pdb_name):
	'''Reformat pdb file so the residue number in column 23-26 is right ordered.
	'''

	command = 'python /home/p/pbryant/pfs/evolution/CATH/parse_pdb_resid.py ' + pdb_name
	out = subprocess.check_output(command, shell = True)#Save parsed pdb
	out = out.decode() #Returns byte
	out = out.split('\n')
	seq = out[0]	
	ca = out[1:-1] 

	reformatted_ca = [] #Save reformatted version

	for line in ca:
		res_number = str(int(line[22:25]))
		new_line = line[0:23]+' '*(3-len(res_number))+res_number+line[26:]
		reformatted_ca.append(new_line)

	#Write to new file
	with open('rf_'+pdb_name, 'w') as file:
		for line in reformatted_ca:
			file.write(line+'\n')

	return None



def write_lddt(indir, out, uid1, uid2):
	'''Writes lddt output to file
	'''
	aln_files = glob.glob(indir +'*.aln')
	for name in aln_files:
		if uid1 in name and uid2 in name:
			pdb_name = name
			break
	#write to file
	lddt_name = pdb_name.split('/')[-1].split('.')[0]+'.lddt'
	with open(lddt_name, 'w') as file:
		out = out.decode() #Returns byte
		file.write(out)

	return None


	
#####MAIN#####
args = parser.parse_args()
indir = args.indir[0]
run_lddt(indir)

