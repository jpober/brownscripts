#!/bin/bash

#SBATCH -J omnical

obsids=("$@")
obs=${obsids[$SLURM_ARRAY_TASK_ID]}
calpar=omni_sol/$obs.pp.fc.npz
python /users/wl42/OmniCal/mp2cal/scripts/omni_run_mwa.py -p 'xx,yy' -C $poscal --ftype='uvfits' --omnipath='./omni_sol/' --wgt_cal --ex_ubls='57_58' --tave $obs
