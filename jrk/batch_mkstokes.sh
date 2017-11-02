#! /bin/bash
#SBATCH -t 0:20:00
#SBATCH -n 2                                                                                                                   
#SBATCH --array=0-999:1%40
###SBATCH --ntasks=1                                                                                                           
#SBATCH --mem=4G                                                                                                              
#SBATCH -J MkStokes
#SBATCH -p jpober-test 
#SBATCH --output=/users/jkerriga/brownscripts/jrk/SlurmOut/FGSub_%A_%a.out 
###SBATCH --output=/users/jkerriga/data/jkerriga/PFHDOutput/fhd_%a/FGSub_%A_%a.out
###SBATCH --error=/users/jkerriga/data/jkerriga/PFHDOutput/fhd_%a/FGSub_%A_%a.err 
source activate PAPER
version=$(($SLURM_ARRAY_TASK_ID + 0))
vsname=''
outdir=/users/jkerriga/data/jkerriga/PSA64FHDLST

obs_list=($(cat obsfits.txt))

echo ${obs_list[$version]}
obs_id=${obs_list[$version]}

echo ${obs_id%.uvfits}[SH]
cd ${outdir}

python ~/brownscripts/jrk/reflag.py  ${obs_id%.uvfits}[SH]
python ~/brownscripts/jrk/mk_stokes.py --stokes='I' ${obs_id%.uvfits}[SH]F



