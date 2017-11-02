#! /bin/bash

#SBATCH -t 0:15:00
#SBATCH -n 1
#SBATCH --array=0-36:1
#SBATCH --mem=2G
###SBATCH -p jpober-test
source activate PAPER
num=$(($SLURM_ARRAY_TASK_ID + 0))
cd ~/data/jkerriga/AllSepsLST/even
declare -a dirt_list
declare -a res_list
dirt_list=(lst*uvH) #Pzen*HPB)
res_list=(lst*uvS) #Pzen*SPB)
#len="${#dirt_list[@]}"
#echo "${dirt_list[@]}"

#for num in $(seq 0 $len);do
#    echo $num
#    echo ${dirt_list[$num]}
python ~/brownscripts/jrk/power_check.py "${dirt_list[$num]}" "${res_list[$num]}"
#done