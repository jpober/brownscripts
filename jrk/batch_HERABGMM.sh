#! /bin/bash
#SBATCH -t 4:00:00
#SBATCH -n 20                                                              
#SBATCH --mem=80G                                                               
#SBATCH -J HERAGaussMix
###SBATCH -p jpober-test
#SBATCH --output=/users/jkerriga/brownscripts/jrk/SlurmOut/WaveletBGMM_%A_%a.out
#SBATCH --array=0-0:1
#outdir=/users/jkerriga/data/jkerriga/PGWideGauss/odd

outdir=/users/jkerriga/data/jkerriga/HERA
cd ${outdir}
source activate PAPER
declare -a dirt_list
declare -a res_list
dirt_list=(zen*UA)
#res_list=(Pzen*SPAB)
lst_list=(lst.*.*uv)
#dirt_list=(lst*uvH)
#res_list=(lst*uvS)


num=$(($SLURM_ARRAY_TASK_ID + 0))
#version=0
#vsname=''

rate=$[($[num%40]*35)/55] # - $[num]%40]
echo $rate
#obs_list=($(cat obsfits.txt))

python ~/brownscripts/jrk/PersistGMM.py #"${dirt_list[@]}" #  "${lst_list[$rate]}"