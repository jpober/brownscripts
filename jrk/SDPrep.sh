#!/bin/bash
#SBATCH -t 1:00:00
#SBATCH -n 2
#SBATCH --array=0-68:1 
#SBATCH --mem=10G                                                                                                                   
#SBATCH -J FHDPrep
#SBATCH --output=SlurmOut/Prep_%A_%a.out
#SBATCH -p jpober-test

source activate PAPER
PSA64Obs='/users/jkerriga/data/jkerriga/PSA64SingleDay'
BrownScripts='/users/jkerriga/brownscripts/jrk'
#cd $PSA64Obs
#obs_list=$(ls -d $PSA64Obs/*O)
obs_list=($(ls -d  $PSA64Obs/zen*O))
echo ${obs_list[*]} |wc -w
cd $BrownScripts
filename=${obs_list[$SLURM_ARRAY_TASK_ID]}
echo $filename
echo 'Rephasing all observations to zenith...'
# Rephase all observations, so that they are phased to zenith
#echo $PSA64Obs/$obs_list
python rephase.py -C psa6240_FHD --onephs $filename

#obs_listM=($(ls -d $PSA64Obs/*M))

echo 'Fixing antenna tables...'
# Fix antenna tables
python fix_anttable.py "${filename}M"

#obs_listT=($(ls -d $PSA64Obs/*T))

echo 'Adding uvw coordinates...'
# Add uvws back onto
python add_uvws.py -C psa6240_FHD "${filename}MT"

#obs_listConvert=($(ls -d $PSA64Obs/*U))
echo 'Converting to UVFITS...'
# Does the miriad to UVFITS conversion while changing NANs to zeros
python miriad2uvfits.py "${filename}MTU"

UVFITS='/users/jkerriga/data/jkerriga/SDOutput'
cd $UVFITS
obs_listFITS=$(ls -d *U.uvfits)
cd $BrownScripts

echo $obs_listFITS > obsfits.txt

#echo 'This is first in the array.'
#obs=( $(cat obs.txt) )
#echo ${obs[100]}