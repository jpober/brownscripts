### Submit the sbatch array command to do omniapply

obs_file_name='/users/wl42/IDL/FHD/Observations/AllNov2016'
#obs_file_name='TESTFC'
mem='60G'
time='2:00:00'
poscal='PhaseII_cal'
#Read the obs file and put into an array, skipping blank lines if they exist
i=0
while read line
do
   if [ ! -z "$line" ]; then
      obs_id_array[$i]=$line
      i=$((i + 1))
   fi
done < "$obs_file_name"

#Create a list of observations using the specified range, or the full observation id file. 
unset good_obs_list
for obs_id in "${obs_id_array[@]}"; do
     good_obs_list+=($obs_id)
done

#Find the number of obsids to run in array
N=${#good_obs_list[@]}                    #Number of files

#N=$(python nfiles.py $paramfile)   	#Generated based on the number of samples per file, the integration time, and the total length

#sbatch -o /dev/null  --array=1-$N --mem=$mem -t $time -n 3 --export=N=$N,fname=$fname,paramfile=$paramfile,calfile=$calfile zeros_job.sh
sbatch -p default-batch --array=0-$(($N-1))%25 --mem=$mem -t $time -n 8 --exclude=node934 --export=N=$N, omniapply.sh ${good_obs_list[@]}
