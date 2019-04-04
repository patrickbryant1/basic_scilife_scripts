#! /usr/bin/env python3
# -*- coding: utf-8 -*-


import argparse
import sys
import matplotlib.pyplot as plt
import numpy as np
import glob
import pandas as pd

#import custom functions
from rnn_input import read_tsv, rmsd_hot
import pdb




#Arguments for argparse module:
parser = argparse.ArgumentParser(description = '''A Recurrent Neural Network for predicting
								RMSD between structural alignments based on sequences from per-residue alignments,
								secondary structure and surface accessibility.''')
 
parser.add_argument('dist_file', nargs=1, type= str,
                  default=sys.stdin, help = 'path to distance file. Format: uid1	uid2	MLdist	RMSD')

parser.add_argument('one_hot_dir', nargs=1, type= str,
                  default=sys.stdin, help = 'path to files with one-hot encodings of alignments, sendary structure and surface acc.')

#Functions

#MAIN
args = parser.parse_args()
dist_file = args.dist_file[0]
one_hot_dir = args.one_hot_dir[0]
#Read tsv
(uids, rmsd_dists) = read_tsv(dist_file, 6)
#Format rmsd_dists into one-hot encoding
rmsd_dists_hot = rmsd_hot(rmsd_dists)


#Get macthing alignments, sendary structure and surface acc
for i in range(0,len(uids)):
	file_name = glob.glob(one_hot_dir + '*/'+uids[i]+'.hot')

	if file_name:
		 np.loadtxt(file_name, dtype=int)


		

	#raise IOerror if file is not found





#When splitting data. Make plots to see rmsd distribution is obtained in train, valid and test