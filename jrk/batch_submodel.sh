#! /bin/bash

###SBATCH -t 0:10:00
###SBATCH -n 1
###SBATCH --array=0-19:1
###SBATCH --mem=2G
###SBATCH -p jpober-test
#source activate PAPER
#num=$(($SLURM_ARRAY_TASK_ID + 0))
#cd ~/data/jkerriga/LastDitchLST
#declare -a dirt_list
#declare -a model_list
#dirt_list=(Pzen*HPA)
#model_list=(Pzen*SPAF)
#len="${#dirt_list[@]}"
#echo "${dirt_list[@]}"

#for num in $(seq 0 $len);do
#    echo $num
#    echo ${dirt_list[$num]}
uv_addsub.py --sub Pzen*HP Pzen*SPF   #"${dirt_list[$num]}" "${model_list[$num]}"