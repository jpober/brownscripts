#! /bin/bash
#SBATCH -t 4:00:00
#SBATCH -n 3                                                            
#SBATCH --mem=4G                                                               
#SBATCH -J HERACNNRFI
###SBATCH -p jpober-test
###SBATCH --output=/users/jkerriga/brownscripts/jrk/SlurmOut/WaveletBGMM_%A_%a.out
#SBATCH --array=0-100:1


outdir=/users/jkerriga/data/jkerriga/HERA2/
cd ${outdir}
source activate PAPER
declare -a dirt_list
declare -a res_list
dirt_list=(zen*uvc)
#res_list=(Pzen*SPAB)
lst_list=(lst.*.*uv)
#dirt_list=(lst*uvH)
#res_list=(lst*uvS)


num=$(($SLURM_ARRAY_TASK_ID + 0))
#version=0
#vsname=''

#rate=$[($[num%40]*35)/55] # - $[num]%40]
#echo $rate
#obs_list=($(cat obsfits.txt))

#python ~/brownscripts/jrk/PersistGMM.py #"${dirt_list[@]}" #  "${lst_list[$rate]}"

#python ~/brownscripts/jrk/AssignPredict.py "${dirt_list[$num]}"
python ~/RFINeural/NeuralAssignPredict.py "${dirt_list[$num]}"