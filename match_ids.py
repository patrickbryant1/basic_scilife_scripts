#! /usr/bin/env python3
# -*- coding: utf-8 -*-


import argparse
import sys
import pdb

#Arguments for argparse module:
parser = argparse.ArgumentParser(description = '''A program that plots evolutionary distance according
	to ML estimations against structural distance in RMSD.''')
 
parser.add_argument('seq_dist_file', nargs=1, type= str,
                  default=sys.stdin, help = 'path to sequence distance file. ')

parser.add_argument('RMSD_file', nargs=1, type= str,
                  default=sys.stdin, help = 'path to file with TMalign output.')



#FUNCTIONS

def read_tsv(tsv_file):
	'''Read tsv file format containing: uid1 \t uid2 \t distance
	'''


	uids = [] #Store uids
	distances = [] #Store distances

	line_number = 0
	with open(tsv_file) as file:
		for line in file:
			if line_number != 0:
				line_number +=1
				line = line.rstrip() #Remove newlines
				line = line.split("\t")
				uid_pair = line[0]+'/'+line[1]
				distance = line[2]

				#For TMalign output
				if len(line) >3:
					chain1 = line[3]
					chain2 = line[4]
					aligned = line[5]
					identity = line[6]

					distance = line[2] + '/' + chain1 + '/' + chain2 + '/' + aligned + '/' + identity
				#Append to lists
				uids.append(uid_pair)
				distances.append(distance)
			else:
				line_number +=1
				continue


	return(uids, distances)



def match_ids(seq_dist_uids, seq_dist_distances, rmsd_uids, rmsd_distances):
	'''A function that matches the ids in each df and fetches the corresponding
	ev_dist and RMSD.
	'''

	sequence_distances = [] #Save sequence distances 
	structural_distances = [] #Save structural distances 

	search_space = list(range(0,len(rmsd_uids))) #Determine intital search space
	
		#Iterate through all seq_dist_uids
	for i in range(0, len(seq_dist_uids)):
		seq_dist_uid_pair = seq_dist_uids[i].split('/') #Split on slash
		seq_dist_uid_1 = seq_dist_uid_pair[0]
		seq_dist_uid_2 = seq_dist_uid_pair[1]

		seq_dist = seq_dist_distances[i]
		sequence_distances.append(seq_dist) #Save seq_dist

		#Match seq_dist_uids to rmsd_uids
		for j in range(0,len(search_space)):
			
			rmsd_uid_pair = rmsd_uids[search_space[j]] #get rmsd_uid_pair for search space index
			if (seq_dist_uid_1 in rmsd_uid_pair) and (seq_dist_uid_2 in rmsd_uid_pair):
				#Save rmsd distance
				structural_distances.append(rmsd_distances[search_space[j]])

				#Print in tsv
				tm_info = rmsd_distances[search_space[j]].split('/')
				print(seq_dist_uid_1 + '\t' + seq_dist_uid_2 + '\t' + seq_dist + '\t' + 
					tm_info[0] + '\t' + tm_info[1] + '\t' + tm_info[2] + '\t' + tm_info[3] + '\t' + tm_info[4])
				search_space.pop(j) #Remove item to reduce search space
				break #When item is found, break out of loop


			else: 
				continue

		


	return(sequence_distances, structural_distances)



#MAIN
args = parser.parse_args()

seq_dist_file = args.seq_dist_file[0]
RMSD_file = args.RMSD_file[0]

#Read tsvs
(seq_dist_uids, seq_dist_distances) = read_tsv(seq_dist_file)
(rmsd_uids, rmsd_distances) = read_tsv(RMSD_file)


#Make sure they have an equal amount of entries
if len(seq_dist_uids) != len(rmsd_uids):
	raise ValueError('The dataframes are not of equal lengths')

#Match uids and
(sequence_distances, structural_distances) = match_ids(seq_dist_uids, seq_dist_distances, rmsd_uids, rmsd_distances)
