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
enc_feature = []
for i in range(0, len(enc1)):
    enc_feature.append(count_aa(enc1[i]))
    enc_feature[i].extend(count_aa(enc2[i]))
#Get RMSDs
rmsds = complete_df['RMSD_x']
bins = np.arange(0,4.5,0.125)

#Bin the TMscore RMSDs
binned_rmsds = np.digitize(rmsds, bins)

#Data
X = np.asarray(enc_feature)
y = np.asarray(binned_rmsds)

#Split
X_train, X_valid, y_train, y_valid = train_test_split(X, y, test_size=0.2, random_state=42)

#RandomForestClassifier
model = RandomForestClassifier(n_estimators=100, bootstrap = True, max_features = 'sqrt')
# Fit on training data
model.fit(X_train, y_train)

#predict
rf_predictions = model.predict(X_valid)
#Average error
average_error = np.average(np.absolute(rf_predictions-y_valid))
print(average_error)
pdb.set_trace()
