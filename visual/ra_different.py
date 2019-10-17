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

parser.add_argument('--calc', nargs=1, type= str,
default=sys.stdin, help = 'either mean or average.')

parser.add_argument('--plot_gradients', nargs=1, type= bool,
default=sys.stdin, help = 'Wether to plot gradients or not.')


def ra_different(df, aln_type, score, cardinalities, calc, plot_num, pdf, fig, ylim, title):
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
            if calc == 'average':
                av= np.average(cut_scores)
            if calc == 'mean':
                av= np.mean(cut_scores)
            avs.append(av)
            js.append(j-step/2)
            total_avs[j-step] = av
            perc_points.append(len(below_df)/len(df)*100)

        plt.plot(js, avs, label = cardinalities[i][1:] +  ' total average', linewidth = 1)
        sizes[cardinality] = [js, perc_points]
    #plt.legend(loc = 'best')
    #plt.xlabel('MLAAdist'+aln_type)
    plt.ylabel(score)
    plt.ylim(ylim)
    plt.xlim([0,6])
    if score == 'RMSD'+aln_type:
        plt.title(title)

    if score == 'DIFF_ACC'+aln_type:
        plot_num+=2
        plt.subplot(plot_num) #set plot_num
        for i in range(len(cardinalities)):
            cardinality = cardinalities[i]
            if cardinality == '_AA20':
                cardinality = ''
            plt.plot(sizes[cardinality][0], sizes[cardinality][1], label = cardinalities[i][1:]+' total average', linewidth = 1)
        plt.xlabel('MLAAdist'+aln_type)
        plt.ylabel('% of points')
        plt.legend(loc = 'best')
        plt.ylim([0,20])
        plt.xlim([0,6])

    return pdf, fig


#####MAIN#####
args = parser.parse_args()
df = pd.read_csv(args.df[0])
outdir = args.outdir[0]
calc = args.calc[0]
plot_gradients = args.plot_gradients[0]

Hgroups = [*Counter(df['H_group']).keys()]
for group in Hgroups:
    partial_df = df[df['H_group']==group]
    i = np.random.randint(len(partial_df), size=1)
    pdb.set_trace()


cardinalities = ['_AA2', '_AA3','_AA6', '_AA20']
pdf = PdfPages(outdir+calc+'_all_cards.pdf')
fig = plt.figure(figsize=(10,10)) #set figsize
plot_nums = [421, 422] #rows,columns,number
aln_types = ['_seqaln', '_straln']
titles = ['Sequence alignments', 'Structural alignments']
ylims = {'RMSD':[0,5], 'DIFFSS':[0, 0.6], 'DIFF_ACC':[0,0.6]}
for i in range(2):
    aln_type = aln_types[i]
    plot_num = plot_nums[i]
    title = titles[i]
    for score in ['RMSD', 'DIFFSS', 'DIFF_ACC']:
        ylim = ylims[score]
        pdf, fig = ra_different(df, aln_type, score, cardinalities, calc, plot_num, pdf, fig, ylim, title)
        plot_num += 2

pdf.savefig(fig)
pdf.close()
