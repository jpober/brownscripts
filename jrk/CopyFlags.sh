#!/bin/bash

cd ~/data/jkerriga/Analysis-10LST


#num=$(($SLURM_ARRAY_TASK_ID + 0))

#dirt_list=($(ls -d ./*HPABR))
#res_list=($(ls -d ./*SPAB))

source activate PAPER
days=('even' 'odd')
for i in ${days[@]};do
    echo $i
    cd ./$i
    dirt_list=($(ls -d ./*HPABR))
    res_list=($(ls -d ./*SPAB))
    for num in $(seq 0 ${#dirt_list[@]});do
	echo $num
	cp -fr ${dirt_list[$num]}/flags ${res_list[$num]}/
	echo ${dirt_list[$num]}/flags ${res_list[$num]}/
    done
    cd ../
done


#cp -r Pzen.2456*[02468].*.uvcRREcACOTUc[HS]PB ./even/
#cp -r Pzen.2456*[13579].*.uvcRREcACOTUc[HS]PB ./odd/