#!/bin/bash
### Submit the sbatch array command to generate a set of miriad files

unset obs_list

basepath=`pwd`

if [ -d "UVFITS" ]; then
	cd UVFITS
	done=($(ls -d *))
	cd ..
fi

if [ ! -d "UVFITS" ]; then
	mkdir UVFITS
fi


if [ $# -eq 0 ];then
	files=($(ls -d *))
else
	files=( "$@" )
fi

for file in ${files[@]}; do
    obsid=`echo $file `
    if [[ " ${done[@]} " =~ " $obsid " ]]; then
    	continue
    else
	echo $obsid
    	obs_list+=($obsid)
    fi
done

nobs=${#obs_list[@]}

sbatch --mem=30G -t 01:10:00 -n 5 --array=0-$(($nobs - 1)) --export=path=$basepath /gpfs_home/alanman/extra_scripts/miriad2uvfits_job.py ${obs_list[@]}

#/gpfs_home/alanman/extra_scripts/miriad_convert.py . ${obs_list[@]}


#sbatch -o /dev/null  --array=1-$N --mem=$mem -t $time -n 3 --export=N=$N,fname=$fname,paramfile=$paramfile,calfile=$calfile zeros_job.sh
#sbatch  --array=1-$N --mem=$mem -t $time -n 3 --export=N=$N,fname=$fname,paramfile=$paramfile,calfile=$calfile zeros_job.sh
