#! /usr/bin/env python3

# -*- coding: utf-8 -*-


import argparse
import sys
import os
import subprocess
#import pexpect
import pdb


#Arguments for argparse module:
parser = argparse.ArgumentParser(description = '''A program that takes aligned residue pairs from the structural
								alignment from TMalign and runs tree-puzzle on them.''')
 
parser.add_argument('align_file', nargs=1, type= str,
                  default=sys.stdin, help = 'path to file with TMalign output.')

#Functions
def get_alignments(file_path):
	'''A function that gets the aligned residue pairs from the structural
	   alignment from TMalign and writes them to files in phylip format
	   including copies of each sequence due to the quartet requirement of tree-puzzle.
	'''

	uid_pairs = [] #List with unique ids
	aligned_seqs = [] #List with aligned residue pairs
	fetch_next_3 = False #Keeping track of lines
	file_names = [] #Store file names

	with open(file_path) as file:
		for line in file:
			if 'Name of Chain' in line:
				line = line.rstrip() #remove \n
				line = line.split("/") #split on /
				uid = line[-1].split(".")[0] #Get uid
				
				uid_pairs.append(uid)



			if fetch_next_3 == True: #Fetch next three lines
				#Fetch lines
				if fetched_lines < 3:
					fetched_lines+=1
					line = line.rstrip() #remove \n
					aligned_seqs.append(line) #Append to list
				else:
					file_name = make_phylip(uid_pairs, aligned_seqs)
					file_names.append(file_name)
					fetch_next_3 = False
					uid_pairs = [] #reset list of pairs
					aligned_seqs = [] #reset aligned seqs
			
				
			if 'denotes aligned residue pairs' in line:
				fetch_next_3 = True
				fetched_lines = 0
			


	return file_names

def make_phylip(uid_pairs, aligned_seqs):
	'''Print phylip format for tree-puzzle calculations
	'''
	#Create text in phylip format
	text = (' 4  ' + str(len(aligned_seqs[0])) + '\n'
			+ uid_pairs[0] + '00|' + aligned_seqs[0] + '\n'
			+ 'copy11111' + '|' + aligned_seqs[0] + '\n'
			+ uid_pairs[1] + '00|' + aligned_seqs[2] + '\n'
			+ 'copy22222' + '|' + aligned_seqs[2] + '\n')
	

	#Define file name
	file_name = uid_pairs[0] + '_' + uid_pairs[1] + '.phy'
	#Open file and write text to it
	with open(file_name, "w") as file:
		file.write(text)

	return file_name


#Main program
args = parser.parse_args()

align_file = args.align_file[0]



#Make .phy files with alignments
file_names = get_alignments(align_file)
print('Number of files to tree-puzzle: '+str(len(file_names)))
#Run tree-puzzle on the files
for name in file_names:
	message_1 = "/home/p/pbryant/pfs/tree-puzzle-5.3.rc16-linux/src/puzzle"
	try:
		p = subprocess.Popen([message_1, name], stdin=subprocess.PIPE)
		p.communicate(b'y\nn\n')[0]
		#p.wait()
	#	p2 = subprocess.Popen(['y'])
	
	except:
		raise IOError(name)
