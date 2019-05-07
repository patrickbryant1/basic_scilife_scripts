#! /usr/bin/env python3
# -*- coding: utf-8 -*-

import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import collections
from scipy import stats
import random

import pdb

#ENCODING INFO

AMINO_ACIDS = { 
'A':1,
'R':2,
'N':3,
'D':4,
'C':5,
'E':6,
'Q':7,
'G':8,
'H':9,
'I':10,
'L':11,
'K':12,
'M':13,
'F':14,
'P':15,
'S':16,
'T':17,
'W':18,
'Y':19,
'V':20,
'X':21, #UNKNOWN
'-':22
}


SECONDARY_STR = {
'G':1,
'H':2,
'I':3,
'T':4,
'E':5,
'B':6,
'S':7,
' ':8,
'-':9
}

#Functions for reading, exploring and visualizing input data

def read_labels(tsv_file):
        '''Read tsv file format containing: uid1 \t uid2 \t ML distance \t RMSD distance
        '''

        distance_dict = {} #Save information for each uid pair
        
        with open(tsv_file) as file:
                for line in file:
                        line = line.rstrip() #Remove newlines
                        line = line.split("\t") #Split on tab
                        uid_pair = (line[0]+'_'+line[1]) #Get uid pair
                        
                        ML_dist = float(line[2]) #aa evolutionary distance from tree-puzzle
                        rmsd_dist = float(line[3]) #rmsd from TMalign
                        Chain1 = int(line[4]) 
                        Chain2 = int(line[5])
                        Aligned = int(line[6])/min(Chain1, Chain2)
                        Identity = float(line[7])
                        

                        distance_dict[uid_pair] = [ML_dist, rmsd_dist, Chain1, Chain2, Aligned, Identity]
                        
                                

        return(distance_dict)


def rmsd_hot(rmsd_dists):
        '''Normalize rmsd distances and convert to one-hot encodings
        '''

        #convert to numpy array
        rmsd_dists = np.array(rmsd_dists)
        #Normalize due to max rmsd
        rmsd_dists = rmsd_dists/max(rmsd_dists)
        #Multiply with 100
        rmsd_dists = rmsd_dists*100
        #Round to 0 decimal pints
        rmsd_dists = np.around(rmsd_dists, decimals = 0)
        #Convert to ints
        rmsd_dists = rmsd_dists.astype(int)
        #Convert to one-hot encoding
        rmsd_dists_hot = np.eye(101)[rmsd_dists]


        return rmsd_dists_hot

def get_locations(encode_locations):
        '''Get all file locations for encodings
        '''

        locations = []
        with open(encode_locations) as file:
                for line in file:
                        line = line.rstrip()
                        locations.append(line)



        return locations


def get_encodings(file_name):
        '''Make a book of all encodings (sentences)
        '''

        aa1 = []
        str1 = []
        acc1 = []
        aa2 = []
        str2 = []
        acc2 = []

        encoding = [aa1, str1, acc1, aa2, str2, acc2] #Save int encoding of sequence

        with open(file_name) as file:
                for line in file:
                        line = line.split(',')
                        aa1.append(AMINO_ACIDS[line[0]])
                        str1.append(SECONDARY_STR[line[1]])
                        acc1.append(int(line[2]))
                        aa2.append(AMINO_ACIDS[line[3]])
                        str2.append(SECONDARY_STR[line[4]])
                        acc2.append(int(line[5]))

                        

                        
        encoding = [aa1, str1, acc1, aa2, str2, acc2]

        return(encoding)

def encoding_distributions(chart_type, encoding_info, title, xlabel, ylabel, bins, out_dir, name, scale, xticks):
        '''Look at distributions of encoded info
        '''
        plt.title(title)
        plt.xlabel(xlabel)
        plt.ylabel(ylabel)

        if chart_type == 'bar':
                counts = collections.Counter(encoding_info)
                #sort dicttionary
                od = collections.OrderedDict(sorted(counts.items()))
                x = []
                y = []
                for key in od:
                        x.append(key)
                        y.append(counts[key])

                plt.bar(x , y, log = scale)

                plt.xticks(x, xticks)


        if chart_type == 'hist':
                plt.hist(encoding_info, bins = bins, log = scale)
                mean = np.sum(encoding_info)/len(encoding_info)
                plt.axvline(mean, color='k', linestyle='dashed', linewidth=1)

        plt.savefig(out_dir+name+'.png')
        plt.close()


def get_labels(encodings, distance_dict, threshold, H_groups):
        '''Get corresponding labels for encodings
        '''

        #Save information in lists for later processing
        uids = [] #save uids
        encoding_list = [] #save encodings
        H_group_list = []
        rmsd_dists = []
        ML_dists = []
        Chains = []
        Align_lens = []
        Identities = []

        letters = []
        structures = []
        accessibilities = []
        seqlens = []
        
        for key in encodings:
                uids.append(key)

                
                [ML_dist, rmsd_dist, Chain1, Chain2, Aligned, Identity] = distance_dict[key]

                if ML_dist <= threshold:
                        #Save all info to see distribtuions and to feed later
                        #From tree-puzzle and TMalign
                        rmsd_dists.append(rmsd_dist)
                        ML_dists.append(ML_dist)
                        Chains.append(Chain1)
                        Chains.append(Chain2)
                        Align_lens.append(Aligned)
                        Identities.append(Identity)

                        #From encoding
                        encoding = encodings[key]


                        for pos in range(0, len(encoding[0])):

                                letters.append(encoding[0][pos])
                                letters.append(encoding[3][pos])
                                structures.append(encoding[1][pos])
                                structures.append(encoding[4][pos])
                                accessibilities.append(encoding[2][pos])
                                accessibilities.append(encoding[5][pos])

                        seqlens.append(len(encoding[0]))
                        #encoding = [np.asarray(enc) for enc in encoding]
                        encoding_list.append(encoding)
                        H_group_list.append(H_groups[key])

        return (uids, encoding_list, rmsd_dists, ML_dists, Chains, Align_lens, Identities, letters, structures, accessibilities, seqlens,  H_group_list)


def word_distributions(encoding_list, bins, out_dir, name):
        '''Count the occurence of all words in all encodings in encoding_list
        '''

        cat_encodings = [j for i in encoding_list for j in i]
                
        plt.hist(cat_encodings, bins = bins, log = True)
        plt.savefig(out_dir+name+'_word_distr.png')
        plt.close()

        return None


def label_distr(x, y, name, out_dir, xlabel, ylabel):
        '''plot distribution of labels (normalized rmsd values)
        '''
        xedges, yedges = np.linspace(0, 9, 10*9), np.linspace(0, 8, 10*8)
        hist, xedges, yedges = np.histogram2d(x, y, (xedges, yedges))

        xidx = np.clip(np.digitize(x, xedges), 0, hist.shape[0]-1)
        yidx = np.clip(np.digitize(y, yedges), 0, hist.shape[1]-1)
        c = hist[xidx, yidx]
        plt.scatter(x, y, c=c)

        #Calculate line of best fit
        (slope, intercept, r_value, p_value, std_err) = stats.linregress(x, y)

        #Desciption
        plt.title('Sequence vs Structural distance' + '\n' + 'R-squared: ' + str((r_value**2).round(3)) +'|' + 'Slope: ' + str(slope.round(3)))
        plt.xlabel(xlabel)
        plt.ylabel(ylabel)
        #Line of best fit
        plt.plot(x, intercept + slope*np.array(x), 'r')
        #Colorbar
        plt.colorbar()

        plt.savefig(out_dir+name+'.png')
        plt.close()
        
        
        return None

def split_on_h_group(encoding_list, H_group_list, unique_groups, counted_groups, percentages, y, out_dir):
	'''Split data so that there are no overlaps btw H-groups in train, validation and test data
	'''
	X_train = []
	y_train = []
	X_valid = []
	y_valid = []
	X_test = []
	y_test = []


	train_groups = []
	

	#1. Shuffle keys of counted_groups, 2 ensures same random shuffling each time
	unique_groups = np.array(unique_groups)
	random.Random(2).shuffle(unique_groups)

	#Count encodings in H-group and assign to split data
	(train_groups, index) = count_encodings(0, encoding_list, unique_groups,counted_groups, percentages[0])

	(valid_groups, index) = count_encodings(index, encoding_list, unique_groups,counted_groups, percentages[1])
	
	test_groups = unique_groups[index:]


	#Assign X and y based on H-group
	for j in range(len(H_group_list)):
		if H_group_list[j] in train_groups:
			#Make check
			if (H_group_list[j] in valid_groups) or (H_group_list[j] in test_groups):
				print('OVERLAP BTW SPLIT ON H-GROUPS!')
				break
			X_train.append(encoding_list[j])
			y_train.append(y[j])
			continue

		if H_group_list[j] in valid_groups:
			#Make check
			if (H_group_list[j] in train_groups) or (H_group_list[j] in test_groups):
				print('OVERLAP BTW SPLIT ON H-GROUPS!')
				break
			X_valid.append(encoding_list[j])
			y_valid.append(y[j])
			continue

		if H_group_list[j] in test_groups:
			#Make check
			if (H_group_list[j] in valid_groups) or (H_group_list[j] in train_groups):
				print('OVERLAP BTW SPLIT ON H-GROUPS!')
				break
			X_test.append(encoding_list[j])
			y_test.append(y[j])
			continue


	
	#Plot RMSD distributions
	encoding_distributions('hist', np.argmax(y_train, axis = 1), 'Distribution of RMSDs for train set. Number of H-groups: ' + str(len(train_groups)), 'Normalized RMSD', 'count', 101, out_dir, 'train_rmsd', False, [])
	encoding_distributions('hist', np.argmax(y_valid, axis = 1), 'Distribution of RMSDs for validation set. Number of H-groups: ' + str(len(valid_groups)), 'Normalized RMSD', 'count', 101, out_dir, 'valid_rmsd', False, [])
	encoding_distributions('hist', np.argmax(y_test, axis = 1), 'Distribution of RMSDs for test set. Number of H-groups: ' + str(len(test_groups)), 'Normalized RMSD', 'count', 101, out_dir, 'test_rmsd', False, [])

	return (X_train, np.asarray(y_train), X_valid, np.asarray(y_valid), X_test, np.asarray(y_test))

	#Add vals of counted groups up to percentage*len(encoding_list)
	#for i in range(len(encoding_list)):

def count_encodings(start, encoding_list, unique_groups, counted_groups, percentage):
	'''Count H-group entries to ensure a suffucient encoding representation
	is present in each split
	'''

	i = start
	count = 0
	save_groups = []
	target = len(encoding_list)*percentage
	while count < target:
		group = unique_groups[i] #
		count += int(counted_groups[group])
		save_groups.append(group)
		i+=1


	return (save_groups, i)


def plot_split(y):
        '''Plots data distributions after train/valid/test split
        '''

        back = np.argmax(y, axis = 1)


        return None


def pad_cut(X, length):
        '''Pad or cut sequences to be of equal lengths and split into 6 different lists
        '''

        X_pad = [[], [], [], [], [], []] #Save padded
        for i in range(len(X)):
    
            if len(X[i][0]) > length:
                for k in range(0, len(X[i])):
                    X_pad[k].append(X[i][k][0:length])
            else:
                pad = [np.pad(inp, (0,length-len(inp)), 'constant') for inp in X[i]]
                
                for l in range(0, len(pad)):
                    X_pad[l].append(pad[l])

        #Convert to arrays
        X_pad = [np.asarray(X_pad[0]), np.asarray(X_pad[1]), np.asarray(X_pad[2]), np.asarray(X_pad[3]), np.asarray(X_pad[4]), np.asarray(X_pad[5])]

        return X_pad
