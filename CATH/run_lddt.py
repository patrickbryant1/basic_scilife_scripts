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

parser.add_argument('outdir', nargs=1, type= str,
default=sys.stdin, help = 'path to output directory.')

parser.add_argument('mode', nargs=1, type= str,
default=sys.stdin, help = 'mode to run lddt in. If mode = "guide", run using aligned pdb files.')

def run_lddt(indir, outdir, mode):
	'''Renumber pdb files to start with atom 0, residue 1.
	'''

	
	#Run LDDT
	if mode == 'guide':
		pdb_files = glob.glob(indir +'*_aln.pdb')
		#Remove files with rf_, these are older versions of the rf
		aln_pdbs = []
		for i in range(len(pdb_files)):
			if 'rf_' not in pdb_files[i]:
				aln_pdbs.append(pdb_files[i])

		pdb_files = aln_pdbs
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
			write_lddt(indir, outdir, out, uid1, uid2)

			#Pop used pdb files
			pdb_files.pop(i)
			pdb_files.pop(0)

	else:
		pdb_files = glob.glob(indir +'*.pdb')
		original_pdb_files = []
		for i in range(len(pdb_files)):
			if '_aln' not in pdb_files[i]:
				original_pdb_files.append(pdb_files[i])
		pdb_files = original_pdb_files
		#Fix pdb files		
		for i in range(len(pdb_files)):
			pdb_name = pdb_files[i].split('/')[-1]
			#Fix pdb files
			move_res_number(pdb_name)
		#Run lddt
		for i in range(len(pdb_files)):
			pdb_name1 = pdb_files[i].split('/')[-1]
			uid1 = pdb_name1.split('.')[0]
			for j in range(1, len(pdb_files)):
				pdb_name2 = pdb_files[j].split('/')[-1]
				
				#Run lddt on fixed pdb files
				try:
					command = 'singularity run --app lDDT /home/p/pbryant/pfs/singularity/ost.img -c -x -t ' + 'rf_'+pdb_name1 + ' rf_'+pdb_name2
					out = subprocess.check_output(command, shell = True)#Save parsed pdb
					uid2 = pdb_name2.split('.')[0]
                	        	#Write to file
					write_lddt(indir, outdir, out, uid1, uid2)
				except:
					print(uid1+'_'+uid2+' failed lddt')
					continue
				
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

	if mode != 'guide':
		res_number = 0
	for line in ca:
		if mode == 'guide':
			res_number = str(int(line[22:26]))
		else:	
			res_number +=1
		new_line = line[0:23]+' '*(3-len(str(res_number)))+str(res_number)+line[26:]
		reformatted_ca.append(new_line)

	#Write to new file
	with open('rf_'+pdb_name, 'w') as file:
		for line in reformatted_ca:
			file.write(line+'\n')

	return None



def write_lddt(indir, outdir, out, uid1, uid2):
	'''Writes lddt output to file
	'''
	aln_files = glob.glob(indir +'*.aln')
	for name in aln_files:
		if uid1 in name and uid2 in name:
			pdb_name = name
			break
	#write to file
	lddt_name = pdb_name.split('/')[-1].split('.')[0]+'.lddt'
	with open(outdir+lddt_name, 'w') as file:
		out = out.decode() #Returns byte
		file.write(out)

	return None


	
#####MAIN#####
args = parser.parse_args()
indir = args.indir[0]
outdir = args.outdir[0]
mode = args.mode[0]
run_lddt(indir, outdir, mode)

