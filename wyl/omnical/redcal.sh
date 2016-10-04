### Submit the sbatch array command to do omnical

obs='1157207048'
poscal='PhaseII_cal_pos'
#pol='xx,yy'
N=1   					#Number of files
mem='80G'
time='10:00:00'

#N=$(python nfiles.py $paramfile)   	#Generated based on the number of samples per file, the integration time, and the total length

#sbatch -o /dev/null  --array=1-$N --mem=$mem -t $time -n 3 --export=N=$N,fname=$fname,paramfile=$paramfile,calfile=$calfile zeros_job.sh
sbatch -p default-batch --array=1-$N --mem=$mem -t $time -n 10 --export=N=$N,obs=$obs,poscal=$poscal, omnical.sh
