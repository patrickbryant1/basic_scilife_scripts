#! /usr/bin/env python3
# -*- coding: utf-8 -*-


import argparse
import sys
import numpy as np
from ast import literal_eval
import pandas as pd

from sklearn.model_selection import train_test_split
from sklearn.metrics import classification_report, confusion_matrix
from sklearn.ensemble import RandomForestClassifier
from collections import Counter

import pdb

#Arguments for argparse module:
parser = argparse.ArgumentParser(description = '''A random RandomForestClassifier for predicting
                                                RMSD between structural alignments.''')

parser.add_argument('dataframe', nargs=1, type= str,
                  default=sys.stdin, help = 'Path to dataframe in .csv.')

parser.add_argument('out_dir', nargs=1, type= str,
                  default=sys.stdin, help = 'Path to output directory. Include /in end')



#MAIN
args = parser.parse_args()
dataframe = args.dataframe[0]
out_dir = args.out_dir[0]

#FUNCTIONS
def count_aa(encoding):
    '''Count differnt amino acid occurence in encoding
    '''
    counts = Counter(encoding)

    aa_feature = []
    for i in range(0,22):
        if i in counts.keys():
            aa_feature.append(counts[i])
        else:
            aa_feature.append(0)

    return aa_feature
#Assign data and labels
#Read df
complete_df = pd.read_csv(dataframe)
#Get MLAAdist
evdist = complete_df['MLAAdist_x']

#Get encodings
enc1 = []
enc2 = []
[enc1.append(literal_eval(x)) for x in complete_df['enc1']]
[enc2.append(literal_eval(x)) for x in complete_df['enc2']]
#Get longest alignment
enc1_feature = []
enc2_feature = []
for i in range(0, len(enc1)):
    enc1_feature.append(count_aa(enc1[i]))
    enc2_feature.append(count_aa(enc2[i]))

pdb.set_trace()
