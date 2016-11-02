#!/bin/bash

#SBATCH -J omnical
i=$SLURM_ARRAY_TASK_ID

python /users/wl42/OmniCal/capo/wyl/scripts/omni_run_multi.py -p 'xx,yy' -C $poscal --ftype='uvfits' --iffits --omnipath='./' $obs
#python omni_run_multi.py -p $pol -C $poscal --ftype='uvfits' --iffits --omnipath='./' $obs
