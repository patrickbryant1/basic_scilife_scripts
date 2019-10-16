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

parser.add_argument('--score', nargs=1, type= str,
default=sys.stdin, help = 'score to plot.')

parser.add_argument('--cardinality', nargs=1, type= str,
default=sys.stdin, help = 'cardinality to plot.')

parser.add_argument('--plot_gradients', nargs=1, type= bool,
default=sys.stdin, help = 'Wether to plot gradients or not.')

parser.add_argument('--plot_percentage', nargs=1, type= bool,
default=sys.stdin, help = 'Wether to plot percentages or not.')

def runnning_average(outdir, complete_df, aln_type, score, cardinality, plot_gradients, plot_percentage, pdf):
    '''Produce running average plots for df
    '''
    figs = plt.figure()
    plot_num = 321 #3 rows 1 columns = 1 graphs per page

    classes = {1.:'Alpha', 2: 'Beta', 3: 'Alpha Beta', 4: 'Few 2ndary structures', 'total': 'Total'}
    colors = {1: 'royalblue', 2: 'k', 3: 'yellow', 4: 'violet', 'total': 'r'}
    ucas = [1.,2.,3.,4.] #unique C.A.s
    fig = plt.figure(figsize=(20,10)) #set figsize
    sizes = {}

    if cardinality == 'AA20':
        cardinality = ''
    score = score+aln_type
    #Plot total average
    avs = [] #Save average score
    js = {} #Save dists
    perc_points = {}
    total_avs = {}
    step = 0.1 #what they used 2009
    df = complete_df
    mldists = np.asarray(df['MLAAdist'+cardinality+aln_type])
    scores = np.asarray(df[score])
    perc_points['total'] = []
    js['total'] = []
    for j in np.arange(min(mldists)+step,max(mldists)+step,step):
        below_df = df[df['MLAAdist'+cardinality+aln_type]<j] #all below j
        below_df = below_df[below_df['MLAAdist'+cardinality+aln_type]>=j-step] #all below j, but at least j-step
        cut_scores = np.asarray(below_df[score])
        av= np.average(cut_scores)
        avs.append(av)
        js['total'].append(j-step/2) #Should be in the middle of interval
        total_avs[j-step/2] = av #used to calculate distance from total average
        perc_points['total'].append(len(below_df)/len(complete_df)*100)

    plt.subplot(plot_num) #set plot_num
    plt.plot(js['total'], avs, label = 'Total average', color = colors['total'], linewidth = 6)

    #Include derivatives
    grads = np.gradient(avs)

    #Save distance from average
    distance_from_av = {}

    for uca in ucas:
        perc_points[uca] = []
        js[uca] = []
        df = complete_df[complete_df['C._x']==uca]
        mldists = np.asarray(df['MLAAdist'+cardinality+aln_type])
        scores = np.asarray(df[score])
        #plt.scatter(mldists, scores, color = 'wheat') #add all points
        avs = [] #Save average score
        distance_from_av[uca] = []

        for j in np.arange(min(mldists)+step,max(mldists)+step,step):
            below_df = df[df['MLAAdist'+cardinality+aln_type]<j]
            below_df = below_df[below_df['MLAAdist'+cardinality+aln_type]>=j-step] #all below j, but at least j-step
            cut_scores = np.asarray(below_df[score])
            if cut_scores.size == 0: #If no scores
                continue
            av= np.average(cut_scores)
            avs.append(av)
            js[uca].append(j-step/2) #Should be in the middle of interval

            #Save av distances and perc points
            #distance_from_av[uca].append(av-total_avs[j-step/2])
            perc_points[uca].append(len(below_df)/len(complete_df)*100)

        perc = np.round(len(df)*100/len(complete_df),2)
        sizes[uca] = perc

        plt.plot(js[uca], avs, label = classes[int(uca)]+', '+str(perc)+' %, '+str(len(df))+' points', color =colors[int(df['C._x'].values[0])], linewidth = 3)


    plt.legend(loc = 'best')
    #plt.xscale('log')
    plt.xlabel('MLAAdist'+cardinality+aln_type)
    plt.ylabel('Running average '+ score)
    plt.title('Running average plot '+aln_type[1:])
    plt.rc('axes', titlesize=18, labelsize=14)
    fig.savefig(outdir+score+cardinality+'.png')#Save

    #Plot gradients
    if plot_gradients == True:
        plot_num+=1
        plt.subplot(plot_num) #set plot_num
        plt.scatter(js['total'][0:40], grads[0:40])
        plt.xlabel('MLAAdist'+cardinality+aln_type)
        plt.ylabel('gradient_'+score)
        plt.title('Running average gradients '+aln_type[1:])
        fig.savefig(outdir+score+cardinality+'_gradients.png')#Save


    if plot_percentage == True:
        plot_num+=1
        plt.subplot(plot_num) #set plot_num
        for label in ['total', 1.,2.,3.,4.]:
            plt.plot(js[label], perc_points[label], label = classes[label], color = colors[label], linewidth = 6)
        plt.xlabel('MLAAdist')
        plt.ylabel('% of points')
        plt.title('Percent of pairs within step')
        plt.legend(loc = 'best')
        #plt.yscale('log')
        fig.savefig(outdir+score+cardinality+'_percentages.png')#Save

    pdf.savefig(fig)
    pdf.close()
    return None


#####MAIN#####
args = parser.parse_args()
df = pd.read_csv(args.df[0])
outdir = args.outdir[0]
aln_type = args.aln_type[0]
score = args.score[0]
cardinality = args.cardinality[0]
plot_gradients = args.plot_gradients[0]
plot_percentage = args.plot_percentage[0]

#Define pdf
pdf = PdfPages(outdir+score+cardinality+'.pdf')
runnning_average(outdir, df, aln_type, score, cardinality, plot_gradients, plot_percentage, pdf)
