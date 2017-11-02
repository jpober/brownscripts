#! /bin/bash
#SBATCH -t 0:45:00
#SBATCH -n 2                                                                                                                                      
#SBATCH --array=0-68:1                                                                                                                                
###SBATCH --ntasks=1                                                                                                                                      
#SBATCH --mem=6G                                                                                                                                        
#SBATCH -p jpober-test 
#SBATCH --output=/users/jkerriga/brownscripts/jrk/SlurmOut/FGSub_%A_%a.out 
###SBATCH --output=/users/jkerriga/data/jkerriga/PFHDOutput/fhd_%a/FGSub_%A_%a.out
###SBATCH --error=/users/jkerriga/data/jkerriga/PFHDOutput/fhd_%a/FGSub_%A_%a.err 
version=$SLURM_ARRAY_TASK_ID
outdir=/users/jkerriga/data/jkerriga/SDOutput
mkdir $outdir/fhd_$version/grid_out

module load ghostscript
module load imagemagick/6.6.4
module load git/2.2.1

ncores=2
FHDpath=$(idl -e 'print,rootdir("fhd")')

echo Using default output directory: $outdir
mkdir -p ${outdir}/fhd_${version}
mkdir -p ${outdir}/fhd_${version}/grid_out
echo Output located at ${outdir}/fhd_${version}

obs_list=($(cat obsfits.txt))

#message=$(sbatch --mem=$mem -t ${wallclock_time} -n ${ncores} --export=ncores=$ncores,outdir=$outdir/fhd_${version},version=$version,thresh=$thresh -o ${outdir}/fhd_${version}/grid_out/firstpass-%A_%a.out -e ${outdir}/fhd_${version}/grid_out/firstpass-%A_%a.err ~/brownscripts/jrk/eor_firstpass_slurm_job.sh ${obs_list})
echo ${obs_list[$version]}
#./eor_firstpass_slurm_job.sh ${obs_list[$version]}
obs_id=${obs_list[$version]}
/usr/local/bin/idl -IDL_DEVICE ps -quiet -IDL_CPU_TPOOL_NTHREADS $ncores -e paper_psa64 -args $obs_id $outdir $version
