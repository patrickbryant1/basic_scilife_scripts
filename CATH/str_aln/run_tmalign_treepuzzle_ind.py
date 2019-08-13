#! /usr/bin/env python3
# -*- coding: utf-8 -*-


import argparse
import sys
import os
import glob
import subprocess
import pdb



#Arguments for argparse module:
parser = argparse.ArgumentParser(description = '''A program that runs TMalign and tree-puzzle and
						receives the resulting output.''')

parser.add_argument('outdir', nargs=1, type= str, default=sys.stdin, help = 'Path to output directory.')
parser.add_argument('hgroup', nargs=1, type= str, default=sys.stdin, help = 'H-group.')
parser.add_argument('puzzle', nargs=1, type= str, default=sys.stdin, help = 'Path to tree-puzzle.')
parser.add_argument('TMalign', nargs=1, type= str, default=sys.stdin, help = 'Path to TMalign.')
parser.add_argument('address', nargs=1, type= str, default=sys.stdin, help = 'Web adress to download from.') #www.cathdb.info/version/v4_2_0/api/rest/id/

#FUNCTIONS
def run_TMalign(outdir, TMalign, hgroup, address):
	'''Download pdb files from CATH api and run TMalign.
	'''
	
	measures = {} #Save RMSD to add with MLAA distance from tree-puzzle
	uids = pd.read_csv(outdir+hgroup+'.txt', sep ='\n')
	#Get pdb files for domains
	for uid in uids:
		#Get pdb file. Make sure they are saved to outdir
		subprocess.call(["wget",address+uid+'.pdb'])			

	#Run TMalign
	for i in range(len(uids)):
		str_i = outdir+uids[i]+'.pdb'
		for j in range(i+1, len(uids)):
			str_j = outdir+uids[i]+'.pdb'
			#Run TMalign and parse output
			tmalign_out = subprocess.check_output([TMalign, str1 , str2]) #Performs optimal structural alignment 
			(tm_aligned_len, rmsd, tmscores, tm_identity, chain_lens, tm_sequences)= parse_tm(tmalign_out)	
			measures[uids[i]+'_'+uids[j]] = [rmsd, tmscores[0], tmscores[1]]
			#Make .phy file of aligned sequences
			make_phylip([uids[i], uids[j]], tm_sequences[0], tm_sequences[1]) 

	return measures, status

def parse_tm(tmalign_out):
	'''A function that gets the uids and the corresponding scores
	and prints them in tsv.
	'''
	
	tmalign_out = tmalign_out.decode("utf-8")
	tmalign_out = tmalign_out.split('\n')
	tmscores = [] #Save TMscores
	for i in range(0, len(tmalign_out)): #Step through all items in list
			
		if 'Aligned length' and 'RMSD' and 'Seq_ID' in tmalign_out[i]:
			row = tmalign_out[i].split(',')
			aligned_len = row[0].split('=')[1].lstrip()
			rmsd = row[1].split('=')[1].lstrip()
			identity = row[2].split('=')[2].lstrip() 
		
		if 'Length of Chain_1:' in tmalign_out[i]:
			len_1 = tmalign_out[i].split(':')[1].split()[0]
				
		if 'Length of Chain_2:' in tmalign_out[i]:
                        len_2 = tmalign_out[i].split(':')[1].split()[0]
		if 'TM-score=' in tmalign_out[i]:
			tmscores.append(tmalign_out[i].split('(')[0].split('=')[1].strip())

	#Get per residue sequence alignments from structural alignment
	sequences = [tmalign_out[-5], tmalign_out[-3]]

	chain_lens = [int(len_1), int(len_2)]
	
			
	return(aligned_len, rmsd, tmscores, identity, chain_lens, sequences)

def make_phylip(uids, query_aln, template_aln):
        '''Print phylip format for tree-puzzle calculations
        '''
        #Create text in phylip format
        text = (' 4  ' + str(len(query_aln)) + '\n'
                        + uids[0] + '00|' + query_aln + '\n'
                        + 'copy11111' + '|' + query_aln + '\n'
                        + uids[1] + '00|' + template_aln + '\n'
                        + 'copy22222' + '|' + template_aln + '\n')


        #Define file name
        file_name = uids[0] + '_' + uids[1] + '.phy'
        #Open file and write text to it
        with open(file_name, "w") as file:
                file.write(text)

        return None

def run_puzzle(indir, puzzle):
        '''Run tree-puzzle and retrieve output
        '''
        for name in glob.glob(indir+"*.phy"): #Use all .phy files
                uid_pairs = name.split('/')[-1].split('.')[0].split('_')
                try:
                        p = subprocess.Popen([puzzle, name], stdin=subprocess.PIPE)
                        p.communicate(b'y\nn\n')[0]
                except:
                        raise IOError(name)


        return None

def parse_puzzle(measures, indir):
	'''Parse output from tree-puzzle and write to dict
	'''
	keys = [*measures] #Make list of keys in dict
	for key in keys:
		uids = key.split('_')
		rmsd, tmscore1, tmscore2 = measures[key] #Get rmsd
		try:
			dist_file = open(indir + key + '.phy.dist')
		except:
			uids = key.split('_')
			dist_file = open(indir + uids[1] + '_' + uids[0] + '.phy.dist')
			measures.pop(key)
			#change key to match other file names
			key = uids[1] + '_' + uids[0]
		for line in dist_file:
			line = line.rstrip()
			line = line.split(" ") #split on double space
			line = list(filter(None, line)) #Filter away empty strings

			if len(line)>2:
				seq_dist = line[-1] #Get ML evolutionary distance between sequences
				measures[key] = [rmsd, tmscore1, tmscore2, seq_dist] 
				break
		dist_file.close()

	return measures


def print_tsv(measures, hgroup):
	'''Print measures in tsv to file
	'''
	with open(hgroup+'.tsv', 'w') as file:
		file.write('uid1\tuid2\tMLAAdist\tRMSD\tTMscore_high\tTMscore_low\n')
		for key in measures:
			uids = key.split('_')
			rmsd, tmscore1, tmscore2, seq_dist = measures[key]
			high_score = max(float(tmscore1), float(tmscore2))
			low_score = min(float(tmscore1), float(tmscore2)) 
			file.write(uids[0]+'\t'+uids[1]+'\t'+seq_dist+'\t'+rmsd+'\t'+str(high_score)+'\t'+str(low_score)+'\n')

	return None

#####MAIN#####
args = parser.parse_args()

outdir = args.outdir[0]
hgroup = args.hgroup[0]
puzzle = args.puzzle[0]
TMalign = args.TMalign[0]
address = args.address[0

#Download pdb files from CATH api and run TMalign
(measures) = run_TMalign(outdir, TMalign)
#Run tree-puzzle on .phy files created from extracted sequence alignments from TMalign struvtural alignments
run_puzzle(outdir, puzzle)
#Parse dist files from tree-puzzle and match to TMalign results
measures = parse_puzzle(measures, outdir)
print_tsv(measures, hgroup)
