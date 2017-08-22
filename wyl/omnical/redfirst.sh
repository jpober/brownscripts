### Submit the sbatch array command to do omnical

obs_file_name='/users/wl42/IDL/FHD/Observations/AllNov2016'
#obs_file_name='ALL'
poscal='PhaseII_cal'
#pol='xx,yy'
mem='60G'
time='1:00:00'

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
#sbatch -o /dev/null  --array=1-$N --mem=$mem -t $time -n 3 --export=N=$N,fname=$fname,paramfile=$paramfile,calfile=$calfile zeros_job.sh
#sbatch -p default-batch --array=0-$(($N - 1)) --mem=$mem -t $time -n 10 --export=N=$N,poscal=$poscal, omnical.sh ${good_obs_list[@]}

sbatch -p default-batch --array=0-$(($N - 1)) --mem=$mem -t $time --exclude=node934 -n 10 --export=N=$N,poscal=$poscal, firstcal.sh ${good_obs_list[@]}


