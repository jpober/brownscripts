#!/bin/bash

#SBATCH -J uvfits_zeros
#SBATCH --mem=50G
#SBATCH -t 8:00:00


#./uvfits_zeros_gen.py hera37_cal.py HERA_sim_2.py -N 94 -o hera37_200ch
./uvfits_zeros_gen.py mwa_128_cal.py HERA_sim_1.py -N 1 -o hera128

