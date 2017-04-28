#!/bin/bash

#SBATCH -J make_noise
i=$SLURM_ARRAY_TASK_ID

outpath=$path"/"$outpath
fname=$path"/"$obs"_"$i


if [ -n "$suffix" ]; then
	/gpfs/data/jpober/alanman/extra_scripts/add_noise.py -o $outpath -s $suffix -N $N $fname
else
	/gpfs/data/jpober/alanman/extra_scripts/add_noise.py -o $outpath -N $N $fname
fi
