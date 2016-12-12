#!/bin/bash

#SBATCH -J omniapply

obsids=("$@")
obs=${obsids[$SLURM_ARRAY_TASK_ID]}

python /users/wl42/OmniCal/capo/wyl/scripts/omni_apply_multi.py -p 'xx,yy' --intype='uvfits' $obs
#python omni_run_multi.py -p $pol -C $poscal --ftype='uvfits' --iffits --omnipath='./' $obs
