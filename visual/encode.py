#! /usr/bin/env python3
# -*- coding: utf-8 -*-

import pdb


#Has to start at 0, since Keras apparently uses the index for lookup
AMINO_ACIDS = { 
'A':0,
'R':1,
'N':2,
'D':3,
'C':4,
'E':5,
'Q':6,
'G':7,
'H':8,
'I':9,
'L':10,
'K':11,
'M':12,
'F':13,
'P':14,
'S':15,
'T':16,
'W':17,
'Y':18,
'V':19,
'X':20, #UNKNOWN
'-':21
}


def get_encodings(file_name):
    '''Make a book of all encodings (sentences)
    Get encoding for file
    '''

    sequences = {} #Save int encoding of sequence
    lses = {} #Save length, start and end
               
    with open(file_name) as file:
        for line in file:
            line = line.rstrip()
            if line[0] == '>':
                uid = line[1:8]
                line = line.split('|')
                lse = line[1]
                lse = lse.split('=')
                l = lse[1].strip('s').strip()
                s = lse[2].strip('e').strip()
                e = lse[3]
                
                lses[uid] = [l,s,e]
                if len(line) > 2:
                    aln_len = line[2].split(':')[1].strip()
                    identity = line[3].split(':')[1].strip()
            else:           
                sequences[uid] = line
                            
    #Get sequences from dict
    keys = [*sequences]
    seq1 = sequences[keys[0]]
    seq2 = sequences[keys[1]]
    
    enc1 = [] #Save encoded sequences
    enc2 = []
    for i in range(0, len(seq1)):
        enc1.append(AMINO_ACIDS[seq1[i]])
        enc2.append(AMINO_ACIDS[seq2[i]])
        
    return(enc1, enc2, int(aln_len), float(identity), lses)