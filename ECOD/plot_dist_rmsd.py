#! /usr/bin/env python3
# -*- coding: utf-8 -*-


import argparse
import sys
import pandas as pd
import subprocess
import pdb
import matplotlib.pyplot as plt
import numpy as np
from scipy import stats
#Arguments for argparse module:
parser = argparse.ArgumentParser(description = '''A program that plots evolutionary distance according
	to ML estimations against structural distance in RMSD.''')
 
parser.add_argument('seq_dist_file', nargs=1, type= str,
                  default=sys.stdin, help = 'path to sequence distance file. ')

parser.add_argument('RMSD_file', nargs=1, type= str,
                  default=sys.stdin, help = 'path to file with TMalign output.')





args = parser.parse_args()

seq_dist_file = args.seq_dist_file[0]
RMSD_file = args.RMSD_file[0]

#Read tsv file as pandas dataframe
seq_dist_df = pd.read_csv(seq_dist_file, sep='\t')#get ev_dist
RMSD_df = pd.read_csv(RMSD_file, sep='\t')#get RMSD
if len(seq_dist_df) != len(RMSD_df):
	raise ValueError('The dataframes are not of equal lengths')


def match_ids(seq_dist_df, RMSD_df):
	'''A function that matches the ids in each df and fetches the corresponding
	ev_dist and RMSD.
	'''

	seq_distances = [] #Save sequence distances 
	structural_distances = [] #Save structural distances 

	end = len(seq_dist_df) #Length of df
	for i in range(0, end):

		uid1 = seq_dist_df['uid1'][i]
		uid2 = seq_dist_df['uid2'][i]
		seq_dist = seq_dist_df['MLdistance'][i]
		seq_distances.append(seq_dist.round(2)) #Add to list and round to tweo decimal points
												#since that is the highest accuracy (accuracy of RMSD from TMalign)


		for j in range(0, end): #Go through df (Exhaustive search, but should not be that big of a problem since the search space is small)
			if RMSD_df['uid1'][j] == uid1 or RMSD_df['uid2'][j] == uid1: #Match to uid1
				if RMSD_df['uid1'][j] == uid2 or RMSD_df['uid2'][j] == uid2: #Match to uid2
					structural_distances.append(RMSD_df['RMSD'][j])


	return(seq_distances, structural_distances)

(seq_distances, structural_distances) = match_ids(seq_dist_df, RMSD_df)

#Calculate stats
(slope, intercept, r_value, p_value, std_err) = stats.linregress(seq_distances, structural_distances)

#Plot
plt.scatter(seq_distances, structural_distances)
plt.title('Defensin related' + '\n' + 'R-squared: ' + str((r_value**2).round(3)))
plt.xlabel('ML AA sequence distance')
plt.ylabel('RMSD')
plt.plot(seq_distances, intercept + slope*np.array(seq_distances), 'r')
plt.show()

#RMSD_df['uid1'] = np.select(RMSD_df['uid2']==ev_dist_df['uid2'])
