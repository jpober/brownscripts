### Submit the sbatch array command to generate a set of uvfits files

fname='ewbase_paper-b.py'
paramfile='EWBASE_sim.py'
calfile='ewbase_paper.py'
N=1   					#Number of files
mem='10G'
time='1:00:00'

N=$(python nfiles.py $paramfile)   	#Generated based on the number of samples per file, the integration time, and the total length

#sbatch -o /dev/null  --array=1-$N --mem=$mem -t $time -n 3 --export=N=$N,fname=$fname,paramfile=$paramfile,calfile=$calfile zeros_job.sh
sbatch -p jpober-test  --array=0-$N --mem=$mem -t $time -n 3 --export=N=$N,fname=$fname,paramfile=$paramfile,calfile=$calfile zeros_job.sh
