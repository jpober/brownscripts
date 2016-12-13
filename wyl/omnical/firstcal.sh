#!/bin/bash

#SBATCH -J firstcal

obsids=("$@")
obs=${obsids[$SLURM_ARRAY_TASK_ID]}

python /users/wl42/OmniCal/capo/wyl/scripts/firstcal_capecod.py -p 'xx,yy' -C $poscal --ftype='uvfits' $obs
#python omni_run_multi.py -p $pol -C $poscal --ftype='uvfits' --iffits --omnipath='./' $obs
