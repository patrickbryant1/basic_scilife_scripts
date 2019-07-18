#! /usr/bin/env python3
# -*- coding: utf-8 -*-


import argparse
import sys
import os
import subprocess
import glob
import pdb

#Arguments for argparse module:
parser = argparse.ArgumentParser(description = '''A program that fixes pdb files
						and then runs lddt on the pdb files.''')

parser.add_argument('indir', nargs=1, type= str,
default=sys.stdin, help = 'path to directory with pdb files.')


def run_lddt(indir, pdb_files):
	'''Run lddt after fixing pdb files
	'''

	
	#Run LDDT
	for i in range(0, len(pdb_files)):
		uid1 = pdb_files[i].split('/')[-1].split('.')[0]
		
		for j in range(i+1, len(pdb_files)):
			uid2 = pdb_files[j].split('/')[-1].split('.')[0]
		
			#Run lddt on fixed pdb files
			command = 'singularity run --app lDDT /home/p/pbryant/pfs/singularity/ost.img -c -x -t ' + 'rf_'+uid1+'.pdb' + ' rf_'+uid2+'.pdb'
			out = subprocess.check_output(command, shell = True)#Save parsed pdb
		
			#Write to file
			write_lddt(indir, out, uid1, uid2)

	return None

def move_res_number(pdb_files):
	'''Reformat pdb file so the residue number in column 23-26 is right ordered.
	'''

        #Fix pdb files
	for pdb_name in pdb_files:
		command = 'python /home/p/pbryant/pfs/evolution/CATH/parse_pdb_resid.py ' + pdb_name
		out = subprocess.check_output(command, shell = True)#Save parsed pdb
		out = out.decode() #Returns byte
		out = out.split('\n')
		ca = out[1:-1] 

		reformatted_ca = [] #Save reformatted version
		new_number=1
		for line in ca:
			res_number = str(int(line[22:26]))
			new_line = line[0:22]+' '*(4-len(str(new_number)))+str(new_number)+line[26:]
			reformatted_ca.append(new_line)
			new_number+=1
		#Write to new file
		outname = 'rf_'+pdb_name.split('/')[-1]
		with open(outname, 'w') as file:
			for line in reformatted_ca:
				file.write(line+'\n')
	
	done = True
	return done



def write_lddt(indir, out, uid1, uid2):
	'''Writes lddt output to file
	'''
	#write to file
	lddt_name = uid1+'_'+uid2+'.lddt'
	with open(lddt_name, 'w') as file:
		out = out.decode() #Returns byte
		file.write(out)
	return None


	
#####MAIN#####
args = parser.parse_args()
indir = args.indir[0]
pdb_files = glob.glob(indir +'*.pdb')
done = move_res_number(pdb_files)
run_lddt(indir, pdb_files)

