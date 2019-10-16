#! /usr/bin/env python3
# -*- coding: utf-8 -*-


import matplotlib.pyplot as plt
from matplotlib.backends.backend_pdf import PdfPages
import pandas as pd
from collections import Counter
import numpy as np
import seaborn as sns
import sys
import argparse

import pdb

#Arguments for argparse module:
parser = argparse.ArgumentParser(description = '''A program that plots running averages.''')

parser.add_argument('--df', nargs=1, type= str,
default=sys.stdin, help = 'path to directory with pdb files.')

parser.add_argument('--outdir', nargs=1, type= str,
default=sys.stdin, help = 'path to output directory.')

parser.add_argument('--aln_type', nargs=1, type= str,
default=sys.stdin, help = '_straln or _seqaln.')

parser.add_argument('--plot_gradients', nargs=1, type= bool,
default=sys.stdin, help = 'Wether to plot gradients or not.')

parser.add_argument('--plot_percentage', nargs=1, type= bool,
default=sys.stdin, help = 'Wether to plot percentages or not.')

def ra_different(df, aln_type, score, cardinalities, plot_num, pdf, fig):
    '''Produce running average plots for df
    '''

    plt.rc('axes', titlesize=10, labelsize=10) #set title and label text sizes
    plt.subplot(plot_num) #set plot_num
    sizes = {} #Save percentage of points in each step
    score = score+aln_type
    for i in range(len(cardinalities)):
        #Plot total average for cardinality
        cardinality = cardinalities[i]
        if cardinality == '_AA20':
            cardinality = ''
        avs = [] #Save average score
        js = [] #Save dists
        perc_points = []
        total_avs = {}
        step = 0.1
        mldists = np.asarray(df['MLAAdist'+cardinality+aln_type])
        scores = np.asarray(df[score])

        for j in np.arange(min(mldists)+step,max(mldists)+step,step):
            below_df = df[df['MLAAdist'+cardinality+aln_type]<j]
            below_df = below_df[below_df['MLAAdist'+cardinality+aln_type]>=j-step]
            cut_scores = np.asarray(below_df[score])
            av= np.average(cut_scores)
            avs.append(av)
            js.append(j-step/2)
            total_avs[j-step] = av
            perc_points.append(len(below_df)/len(df)*100)

        plt.plot(js, avs, label = cardinalities[i][1:] +  ' total average', linewidth = 1)
        sizes[cardinality] = [js, perc_points]
    plt.legend(loc = 'best')
    plt.xlabel('MLAAdist'+aln_type)
    plt.ylabel('Running average '+ score)
    plt.title('Running average plot '+score+' '+aln_type[1:])

    if score == 'DIFF_ACC':
        plt.subplot(plot_num+1) #set plot_num
        for i in range(len(cardinalities)):
            cardinality = cardinalities[i]
            if cardinality == '_AA20':
                cardinality = ''
            plt.plot(sizes[cardinality][0], sizes[cardinality][1], label = cardinalities[i][1:]+' total average', linewidth = 1)
        plt.xlabel('MLAAdist'+aln_type)
        plt.ylabel('% of points')
        plt.legend(loc = 'best')

    return pdf, fig


#####MAIN#####
args = parser.parse_args()
df = pd.read_csv(args.df[0])
outdir = args.outdir[0]
aln_type = args.aln_type[0]
plot_gradients = args.plot_gradients[0]
plot_percentage = args.plot_percentage[0]

cardinalities = ['_AA2', '_AA3','_AA6', '_AA20']
pdf = PdfPages(outdir+aln_type[1:]+'.pdf')
fig = plt.figure(figsize=(20,10)) #set figsize
plot_num = 411 #rows,columns,number

for score in ['RMSD', 'DIFFSS', 'DIFF_ACC']:
    pdf, fig = ra_different(df, aln_type, score, cardinalities, plot_num, pdf, fig)
    plot_num += 1

pdf.savefig(fig)
pdf.close()
