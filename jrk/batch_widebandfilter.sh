#! /bin/bash
#SBATCH -t 0:30:00
#SBATCH -n 2                                                                                                            
#SBATCH --mem=10G                                                                                                                                             
#SBATCH -J WidebandFilter                                                                                                                                     
#SBATCH -p jpober-test                                                                                                                                    
#SBATCH --output=/users/jkerriga/brownscripts/jrk/SlurmOut/WideBand_%A_%a.out 
#SBATCH --array=0-999:1
outdir=/users/jkerriga/data/jkerriga/AnalysisLST

cd ${outdir}
source activate PAPER
declare -a dirt_list
declare -a res_list
dirt_list=(Pzen*HPA)
res_list=(Pzen*SPA)

#dirt_list=(lst*uvHB)
#res_list=(lst*uvSB)


num=$(($SLURM_ARRAY_TASK_ID + 3999))
#version=0
#vsname=''


#obs_list=($(cat obsfits.txt))

#echo ${obs_list[$version]}
#obs_id=${obs_list[$version]}

#echo ${obs_id%.uvfits}[SH]P
#cd ${outdir}

python ~/capo/pspec_pipeline/pspec_prep.py -C psa6240_FHD -a cross --nogain --nophs --clean=1e-9 --horizon=15 --window='blackman-harris' "${dirt_list[$num]}" "${res_list[$num]}"

python ~/capo/dcj/scripts/xrfi_simple.py -n 3 "${dirt_list[$num]}""B" "${res_list[$num]}""B"

#python ~/brownscripts/jrk/FilterHorizon.py "${dirt_list[$num]}" "${res_list[$num]}"

#Pzen.2456*${SLURM_ARRAY_TASK_ID}.*.uvcRREcACOTUc[HS]P #${obs_id%.uvfits}[SH]P

#python ~/capo/zsa/scripts/pspec_prep.py -C psa6240_FHD -a cross --nogain --nophs --clean=1e-9 --horizon=15 --window='blackman-harris' Pzen.2456*${SLURM_ARRAY_TASK_ID}.*.uvcRREcACOTUcSPA

#python ~/capo/zsa/scripts/pspec_prep.py -C psa6240_FHD -a cross --window=blackman-harris --nogain --nophs --clean=1e-9 --horizon=15 Pzen.2456*${SLURM_ARRAY_TASK_ID}.*.uvcRREcACOTUcHPA
