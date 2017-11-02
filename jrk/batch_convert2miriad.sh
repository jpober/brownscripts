#! /bin/bash
#SBATCH -t 1:00:00
#SBATCH -n 2                                                                                                                                      
#SBATCH --array=0-88:1                                                                                                                             
###SBATCH --ntasks=1                                                                                                                  
#SBATCH --mem=4G                                                                                                                                  
#SBATCH -p jpober-test 
#SBATCH --output=/users/jkerriga/brownscripts/jrk/SlurmOut/FHD2Miriad_%A_%a.out 
version=$(($SLURM_ARRAY_TASK_ID + 0))
outdir=/users/jkerriga/data/jkerriga/2DayOutput
#mkdir $outdir/fhd_$version/grid_out
flag='-s'
module load ghostscript
module load imagemagick/6.6.4
module load git/2.2.1
echo Output located at ${outdir}/fhd_${version}

cd ${outdir}/fhd_${version}/vis_data
python ~/brownscripts/jrk/sav2miriad.py Pzen.* ../metadata/Pzen.* ${flag}

if [ $flag == "-s" ]; then
    cp -rf *S /users/jkerriga/data/jkerriga/FGLST/
    rm -r *S
else
    cp -rf *H /users/jkerriga/data/jkerriga/FGLST/
    rm -r *H
fi
