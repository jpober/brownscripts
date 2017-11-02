#! /bin/bash
#SBATCH -t 0:20:00
#SBATCH -n 2                                                                   
#SBATCH --mem=40G                                                                  
#SBATCH -J WidebandFilter                                                                              
###SBATCH -p jpober-test
                                                                                                      
                             
#SBATCH --output=/users/jkerriga/brownscripts/jrk/SlurmOut/WideBand_%A_%a.out 
#SBATCH --array=0-999:1
#outdir=/users/jkerriga/data/jkerriga/PGWideGauss/odd
outdir=/users/jkerriga/data/jkerriga/PGBH/even
cd ${outdir}
source activate PAPER
declare -a dirt_list
declare -a res_list
dirt_list=(Pzen*HPAB)

lst_list=(lst.*.*uvHB)
#dirt_list=(lst*uvH)
#res_list=(lst*uvS)


num=$(($SLURM_ARRAY_TASK_ID + 1999))
#version=0
#vsname=''

rate=$[($[num%40]*35)/40] # - $[num]%40]
echo $rate
#obs_list=($(cat obsfits.txt))

#echo ${obs_list[$version]}
#obs_id=${obs_list[$version]}

#echo ${obs_id%.uvfits}[SH]P
#cd ${outdir}

#python ~/brownscripts/jrk/pspec_prep.py -C psa6240_FHD -a cross --nogain --nophs --clean=1e-9 --horizon=15 --window='blackman-harris' "${dirt_list[$num]}" "${res_list[$num]}"

#python ~/capo/dcj/scripts/xrfi_simple.py -n 3 "${dirt_list[$num]}""B" "${res_list[$num]}""B"

#python ~/brownscripts/jrk/FilterHorizon.py "${dirt_list[$num]}" "${res_list[$num]}"

#Pzen.2456*${SLURM_ARRAY_TASK_ID}.*.uvcRREcACOTUc[HS]P #${obs_id%.uvfits}[SH]P

#python ~/capo/zsa/scripts/pspec_prep.py -C psa6240_FHD -a cross --nogain --nophs --clean=1e-9 --horizon=15 --window='gaussian0.4' Pzen.2456*${SLURM_ARRAY_TASK_ID}.*.uvcRREcACOTUcSPA

#python ~/capo/pspec_pipeline/pspec_prep.py -C psa6240_FHD -a cross --window='blackman-harris' --nogain --nophs --clean=1e-9 --horizon=15 Pzen.2456*${SLURM_ARRAY_TASK_ID}.*.uvcRREcACOTUc[HS]PA

#python ~/capo/pspec_pipeline/pspec_prep.py -C psa6240_FHD -a cross --window='blackman-harris' --nogain --nophs --clean=1e-8 --horizon=15 Pzen.2456*${SLURM_ARRAY_TASK_ID}.*.uvcRREcACOTUcSPA

python ~/brownscripts/jrk/WaveletBGMM.py "${dirt_list[$num]}"  "${lst_list[$rate]}"