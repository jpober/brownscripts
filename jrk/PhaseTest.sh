#!/bin/bash
#SBATCH -t 5:00:00                                                                                                                                          
#SBATCH -n 1
#SBATCH --array=3-7:1 
#SBATCH --mem=6G                                                                                                                                            
#SBATCH -J PhsTest
#SBATCH --output=SlurmOut/PHSTest_%A_%a.out
####SBATCH -p jpober-test

source activate PAPER
PSA64Obs='/users/jkerriga/data/jkerriga/PSA64TwoDay'
outpath='/users/jkerriga/data/jkerriga/FGLST'
BrownScripts='/users/jkerriga/brownscripts/jrk'
#cd $PSA64Obs
#obs_list=$(ls -d $PSA64Obs/*O)

cd $outpath
time=$SLURM_ARRAY_TASK_ID
echo $time
echo 'Rephasing all observations to zenith...'
python ${BrownScripts}/phasetest.py ${PSA64Obs}/zen.*.${SLURM_ARRAY_TASK_ID}*.uvcRREcACOTUcP



