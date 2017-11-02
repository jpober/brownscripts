#! /bin/bash
#SBATCH -t 2:00:00
#SBATCH -n 12                                                                                                                                     
#SBATCH --array=0-0:1%25
###SBATCH --ntasks=1                                                                                                                                      
#SBATCH --mem=60G                                                                                                                                   
#SBATCH -p jpober-test 
##SBATCH --output=/users/jkerriga/brownscripts/jrk/SlurmOut/FGSub_%A_%a.out 
###SBATCH --output=/users/jkerriga/data/jkerriga/PFHDOutput/fhd_%a/FGSub_%A_%a.out
###SBATCH --error=/users/jkerriga/data/jkerriga/PFHDOutput/fhd_%a/FGSub_%A_%a.err 
source activate PAPER
#version=$(($SLURM_ARRAY_TASK_ID + 0))
version='Deconvolve'
vsname=''
outdir=/users/jkerriga/data/jkerriga/2DayOutput

module load ghostscript
module load imagemagick/6.6.4
module load git/2.2.1

ncores=2
FHDpath=$(idl -e 'print,rootdir("fhd")')

echo Using default output directory: $outdir
mkdir -p ${outdir}/fhd_${vsname}${version}
mkdir -p ${outdir}/fhd_${vsname}${version}/grid_out
echo Output located at ${outdir}/fhd_${vsname}${version}

#obs_list=($(cat Analysis_obsfits.txt))


echo ${obs_list[$version]}
#obs_id=${obs_list[$version]}
obs_id='Pzen.2456242.40348.uvcRREcACOTUc.uvfits'
/usr/local/bin/idl -IDL_DEVICE ps -quiet -IDL_CPU_TPOOL_NTHREADS $ncores -e test_psa64 -args $obs_id $outdir ${vsname}${version}

#cd ${outdir}/fhd_${vsname}${version}/vis_data
#python ~/brownscripts/jrk/sav2miriad.py Pzen.* ../metadata/Pzen.* -s
#cp -rf *UcS ../../../PSA64FHDLST/
#cp -rf *UcH ../../../PSA64FHDLST/
#rm -r *UcS
#rm -r *UcH

