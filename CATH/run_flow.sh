#!/bin/bash -l


#A bash script for getting CATH files for a certain H-group.
#Sequence and structural alignments are then performed.
#Different metrics are thereafter calculated: secondary structure,
#surface accessibility, RMSD, Maximum likelihood amino acid distance,
#TMscore, LDDT-score

#Filter the clustered sequences according to experimental method used for
#structural determination and X-ray resolution
/home/pbryant/evolution/CATH/pdb_filter.py /home/pbryant/data/CATH/mmseqs2/cath-domain-seqs-0.2_rep_seq.fasta /home/pbryant/data/CATH/Xray_res2.6Å.txt /home/pbryant/data/CATH/mmseqs2/



$FILE_NAME="H-group"

#Make results directory
mkdir $RESULTS_DIR
#Go to results directory
cd $RESULTS_DIR
#Program input paths
CATH_API=www.cathdb.info/version/v4_2_0/api/rest/id/
HHBLITS=/home/p/pbryant/pfs/hh-suite/build/bin/hhblits
HHALIGN=/home/p/pbryant/pfs/hh-suite/build/bin/hhalign
UNIPROT=/home/p/pbryant/pfs/uniprot20_2016_02/data/uniprot20_2016_02

PUZZLE=/home/p/pbryant/pfs/tree-puzzle-5.3.rc16-linux/src/puzzle
TMALIGN=/home/p/pbryant/pfs/TMalign
TMSCORE=/home/p/pbryant/pfs/TMscore
#Get ids and run hhalign
/home/p/pbryant/pfs/evolution/CATH/get_data.py /home/p/pbryant/pfs/data/CATH/mmseqs2/h_grouped/below95_2.6Å_above2_grouped/$FILE_NAME/ $RESULTS_DIR/ $HHBLITS $HHALIGN $UNIPROT 15 $CATH_API

#Run dssp
DSSP=/home/p/pbryant/pfs/dssp
wait
mkdir $RESULTS_DIR/dssp
wait
#Run dssp for pdbs
for file in $RESULTS_DIR/*.pdb
do
if [[ $file != *aln* ]] && [[ $file != *rf_* ]];  #If no aln or rf in name
then
$DSSP $file > $file'.dssp'
fi
done
wait
#Move all files to dssp dir
mv $RESULTS_DIR/*.dssp $RESULTS_DIR/dssp/

wait
#Make TMalign dir
mkdir $RESULTS_DIR/TMalign/
wait
#Run TMalign and tree-puzzle
singularity exec /pfs/nobackup/home/p/pbryant/singularity/bio.sif /pfs/nobackup/home/p/pbryant/evolution/CATH/run_tmalign_treepuzzle_ind.py $RESULTS_DIR/ $RESULTS_DIR/TMalign/ /home/p/pbryant/pfs/data/CATH/below95_2.6Å_above2_grouped/$FILE_NAME/ $FILE_NAME $PUZZLE $TMALIGN

#Go to where the files are
cd $RESULTS_DIR/TMalign/
wait
#Run lddt for pdbs
/home/p/pbryant/pfs/evolution/CATH/run_lddt.py $RESULTS_DIR/TMalign/ $RESULTS_DIR/TMalign/ guide

wait
#Make TMScore dir
mkdir $RESULTS_DIR/TMscore/
cd $RESULTS_DIR
wait
#Run TMscore and tree-puzzle
singularity exec /pfs/nobackup/home/p/pbryant/singularity/bio.sif /pfs/nobackup/home/p/pbryant/evolution/CATH/run_tmscore_treepuzzle.py $RESULTS_DIR/ $RESULTS_DIR/TMscore/ /home/p/pbryant/pfs/data/CATH/below95_2.6Å_above2_grouped/$FILE_NAME/ $FILE_NAME $PUZZLE $TMSCORE

wait
#Run lddt for aligned pdbs
/home/p/pbryant/pfs/evolution/CATH/run_lddt.py $RESULTS_DIR/ $RESULTS_DIR/TMscore/ guide
