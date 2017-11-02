#! /bin/bash
#SBATCH -t 3:00:00
#SBATCH -n 1                                                                  
                                      
#SBATCH --mem=10G                                                               
                                                                              
#SBATCH -J WidebandFilter                                                       
                                                                              
###SBATCH -p jpober-test                                                       
#SBATCH --output=/users/jkerriga/brownscripts/jrk/SlurmOut/WideBand_%A_%a.out 
#SBATCH --array=243-365:2
outdir=/users/jkerriga/data/jkerriga/PaperData/oddIH
modeldir=/users/jkerriga/data/jkerriga/PGModels
#cd ${outdir}
#SLURM_ARRAY_TASK_ID=242
source activate PAPER
declare -a dirty_list
declare -a res_list
declare -a model_list

##Get models and filter
#cd $modeldir
#model_list=(Pzen.2456${SLURM_ARRAY_TASK_ID}.*.*MPA)
#python ~/brownscripts/jrk/model_prep.py -C psa6240_FHD -a cross --nogain --nophs --clean=1e-10 --horizon=0 --model --window='blackman-harris' "${model_list[@]}"
#python ~/brownscripts/jrk/ModelDelayFilter.py Pzen.2456${SLURM_ARRAY_TASK_ID}.*.*MPAF

##Move to work directory and subtract models
cd ${outdir}
dirty_list=(Pzen.2456${SLURM_ARRAY_TASK_ID}.*.*HPA)
#modelF_list=(/users/jkerriga/data/jkerriga/PGModels/Pzen.2456${SLURM_ARRAY_TASK_ID}.*.*MPAFF)
#python ~/brownscripts/jrk/uv_sub.py --dirty="Pzen.2456${SLURM_ARRAY_TASK_ID}.*.*HPA" --model="/users/jkerriga/data/jkerriga/PGModels/Pzen.2456${SLURM_ARRAY_TASK_ID}.*.*MPAFF" -s 'SPA'

##Wideband Filter
res_list=(Pzen.2456${SLURM_ARRAY_TASK_ID}.*.*SPA)
#python ~/brownscripts/jrk/pspec_prep.py -C psa6240_FHD -a cross --nogain --nophs --clean=1e-9 --horizon=15 --window='blackman-harris' "${dirty_list[@]}"

python ~/capo/pspec_pipeline/pspec_prep.py -C psa6240_FHD -a cross --nogain --nophs --clean=1e-9 --horizon=-15  --window='blackman-harris' "${dirty_list[@]}" "${res_list[@]}"

#python ~/brownscripts/jrk/PostWideBandFilter.py Pzen.2456${SLURM_ARRAY_TASK_ID}.*.*HPAB Pzen.2456${SLURM_ARRAY_TASK_ID}.*.*SPAB

python ~/brownscripts/jrk/xrfi_simple_jrk.py -n 3 Pzen.2456${SLURM_ARRAY_TASK_ID}.*.*HPAB Pzen.2456${SLURM_ARRAY_TASK_ID}.*.*SPAB #"${dirty_list[@]}""B" "${res_list[@]}""B"

