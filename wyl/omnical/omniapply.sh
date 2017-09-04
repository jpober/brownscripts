#!/bin/bash

#SBATCH -J omniapply

obsids=("$@")
obs=${obsids[$SLURM_ARRAY_TASK_ID]}

python /users/wl42/OmniCal/mp2cal/scripts/omni_apply_mwa.py -p 'xx,yy' -C 'PhaseII_cal' --omnipath='./omni_sol/' --metafits='/users/wl42/data/wl42/Nov2016EoR0/' --ave --intype='uvfits' --appfhd $obs

