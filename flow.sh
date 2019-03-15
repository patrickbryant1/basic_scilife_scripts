#!/usr/bin/env bash

#A script for running a complete round of pairwise structural alignment with TMalign of 
#.pdb structure files followed by RMSD score parsing and residue pair alignment extraction.
#The extracted residue alignments are then used for calculateng sequence distance
#using Tree-puzzle.


#Run TM-align to get RMSDs
/home/pbryant/basic_scilife_scripts/get_rmsd.py /home/pbryant/data/ECOD/defensin_related/defensin_related_pdb_id.txt /home/pbryant/data/pdb/defensin_related/ > defensin_rmsd.txt
wait

#Parse oputput into .tsv with uids and RMSD between uid1 and uid2. Written to file defensin_rmsd.txt
/home/pbryant/basic_scilife_scripts/parse_rmsd.py /home/pbryant/data/ECOD/defensin_related/defensin_related_pdb_id.txt defensin_rmsd.txt > defensin_rmsd.tsv
wait

#Get the alignments from TMalign and run tree-puzzle on them.
#Returns files in the form of uid1_uid2.[phy].[dist, puzzle]
/home/pbryant/basic_scilife_scripts/straln_to_phylip.py /home/pbryant/data/ECOD/defensin_related/defensin_related_pdb_id.txt defensin_rmsd.txt
wait
#Get all pairwise sequence distances from tree-puzzle into tsv
/home/pbryant/basic_scilife_scripts/parse_dist.py /home/pbryant/results/defensin_related/complete_workflow/20190315/ > defensin_related.dist.tsv
wait

#Plot ML sequence distance from tree-puzzle against RMSD from TMalign
/home/pbryant/basic_scilife_scripts/plot_dist_rmsd.py ./defensin_related.dist.tsv ./defensin_rmsd.tsv
