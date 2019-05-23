#! /usr/bin/env python3
# -*- coding: utf-8 -*-


import argparse
import sys
import os
import glob
import subprocess
import pdb



#Arguments for argparse module:
parser = argparse.ArgumentParser(description = '''A program that runs TMalign and tree-puzzle and
						receives the resulting output.''')
 
parser.add_argument('puzzle', nargs=1, type= str, default=sys.stdin, help = 'Path to tree-puzzle.')

parser.add_argument('TMalign', nargs=1, type= str, default=sys.stdin, help = 'Path to TMalign.')


#FUNCTIONS


def run_puzzle(puzzle):
	'''Run tree-puzzle and retrieve output
	'''
	for name in glob.glob("*.phy"): #Use all .phy files
		uid_pairs = name.split('/')[-1].split('.')[0].split('_')
		try:
			p = subprocess.Popen([puzzle, name], stdin=subprocess.PIPE)
			p.communicate(b'y\nn\n')[0]
			output = p.stdout.readline()
			pdb.set_trace()
		except:
			raise IOError(name)


	return None

#####MAIN#####
args = parser.parse_args()

puzzle = args.puzzle[0]
TMalign = args.TMalign[0]

run_puzzle(puzzle)
