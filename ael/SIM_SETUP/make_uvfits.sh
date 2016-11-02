### Submit the sbatch array command to generate a set of uvfits files

fname='HERA37'
paramfile='HERA_sim_3.py'
calfile='hera37_cal.py'
N=1   					#Number of files
mem='10G'
time='00:30:00'

N=$(python nfiles.py $paramfile)   	#Generated based on the number of samples per file, the integration time, and the total length

#sbatch -p jpober-test  -o "slurm-"$fname"_%a.out" --array=0-$N --mem=$mem -t $time -n 1 --export=N=$N,fname=$fname,paramfile=$paramfile,calfile=$calfile zeros_job.sh
sbatch  --array=0-$N -o "slurm-"$fname"_%a.out" --mem=$mem -t $time -n 1 --export=N=$N,fname=$fname,paramfile=$paramfile,calfile=$calfile zeros_job.sh

