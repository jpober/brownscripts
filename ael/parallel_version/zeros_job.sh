#!/bin/bash

#SBATCH -J uvfits_zeros
i=$SLURM_ARRAY_TASK_ID

python uvfits_zeros_par.py $calfile $paramfile -i $i -N $N -o `pwd `'/'$fname
