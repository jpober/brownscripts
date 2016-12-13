#!/bin/bash

#SBATCH -J omnical

obsids=("$@")
obs=${obsids[$SLURM_ARRAY_TASK_ID]}
calpar=$obs.pp.fc.npz
python /users/wl42/OmniCal/capo/wyl/scripts/omni_run_multi.py -p 'xx,yy' -C $poscal --ftype='uvfits' --iffits --omnipath='./omni_sol/' --calpar=$calpar $obs
#python omni_run_multi.py -p $pol -C $poscal --ftype='uvfits' --iffits --omnipath='./' $obs
