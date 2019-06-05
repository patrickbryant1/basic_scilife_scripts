#! /usr/bin/env python3
# -*- coding: utf-8 -*-


import argparse
import sys
import matplotlib.pyplot as plt
import numpy as np
from scipy import stats


import pdb




#Arguments for argparse module:
parser = argparse.ArgumentParser(description = '''A program that plots evolutionary distance according
	to ML estimations against structural distance in RMSD.''')
 
#parser.add_argument('dist_file', nargs=1, type= str,
#                  default=sys.stdin, help = 'path to distance file. Format: uid1	uid2	MLdist	RMSD')



#FUNCTIONS

def read_tsv(tsv_file, threshold):
	'''Read tsv file format containing: uid1 \t uid2 \t ML distance \t RMSD distance
	'''


	ML_dists = [] #Store ML distances
	rmsd_dists = [] #Store rmsd distances
	Z = [] #Store x,y vals for surface plot.

	with open(tsv_file) as file:
		for line in file:
            print(line)
            pdb.set_trace()
            if 'uid1' not in line:
                ML_dist = round(float(line[2]), 2)
                if ML_dist <= threshold:
                    ML_dists.append(ML_dist)
                    rmsd_dists.append(float(line[3]))
                    Z.append((ML_dist,float(line[3])))
                else:
                    continue
                
            else:
                continue

	
	return(ML_dists, rmsd_dists, Z)

##def scatter_plot(X, Y):
##	'''
##	'''
##
##	(slope, intercept, r_value, p_value, std_err) = stats.linregress(X, Y)
##
##
##	#Plot
##	plt.scatter(X, Y)
##	plt.title('Sequence vs Structural distance' + '\n' + 'R-squared: ' + str((r_value**2).round(3)))
##	plt.xlabel('ML AA sequence distance')
##	plt.ylabel('RMSD')
##	plt.plot(X, intercept + slope*np.array(X), 'r')
##	plt.show()
##
##	return None


#MAIN
#args = parser.parse_args()
#dist_file = args.dist_file[0]

#Read tsv
#(ML_dists, rmsd_dists, Z) = read_tsv(dist_file, 6)



