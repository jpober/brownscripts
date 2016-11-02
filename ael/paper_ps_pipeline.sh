#!/bin/bash

## A wrapper to run the PAPER power spectrum pipeline on FHD simulated data.

#Parse input arguments. So far, only option is "c" for configuration file. All other arguments passed are the MIRIAD files to process.
touch tmp
skiptoplot=0
while [[ $@ ]]
do
	OPTIND=1
	if [[ $1 =~ ^- ]]
	then
	    getopts ":c:s" option
	    case $option in
		c) configfile="$OPTARG"
		   shift ;;
		s) skiptoplot=1
	    esac
	else
	    echo $1 >> tmp
	fi
	shift
done


fileItemString=$(cat tmp |tr "\n" " ")
infiles=( $fileItemString )
rm tmp

. $configfile

#Ensure the calfile is nearby
calfile=( $(find . -name $cal".py" -type 'f') )
if [ ${#calfile[@]} -eq 0 ]; then
	calfile=( $(find -L ~/SIM_SETUP -name $cal".py" -type 'f') )    #-L is needed to follow symbolic links
	if [ ! ${#calfile[@]} -eq 0 ]; then
		cp ${calfile[0]} .
	else
		echo "Error: Cal file not found. Exiting."
		exit
	fi
fi
   
if [ $skiptoplot -eq 0 ]; then

   #calfile='psa6622_v003_paper128'
   
   dir=$(dirname ${infiles[0]})
   olddir=$(pwd)
   cd $dir
   
   
   files=()
   for f in ${infiles[@]};
   do
   #	count=$(python -c "print \""$f"\".count('.')")
   #	if [ $count -eq 0 ]; then
   	if [[ ! $f == *'.'* ]]; then
   		cp -r $f $f'.'   #Important --> Needs a suffix
   		files=( ${files[@]} $f'.' )
   	else
   		files=( ${files[@]} $f )
   	fi
   done
   
   reducfiles=()
   for f in ${files[@]};
   do
   	f=$(basename $f)
   	tmp=( $(find . -name $f'A' -type 'd') )
   	if [ ${#tmp[@]} -eq 0 ]; then
   		reducfiles=($f ${reducfiles[@]})
   	fi
   done
   
   if [ ! ${#reducfiles[@]} -eq 0 ]; then
   	echo 'Reduce_Seps'
   	Reduce_Seps128.sh ${reducfiles[@]}
   #	files=$(find . -name "*A" -type 'd' -mtime -10)   #Find MIRIAD files made less than 10 seconds ago ending in A
   fi
   
   #See if stokes-I files exist. If not, run mk_stokes
   stokefiles=()
   for f in ${files[@]};
   do
   	f=$(basename $f)
   	tmp=( $(find . -name $f'AP' -type 'd') )
   	if [ ${#tmp[@]} -eq 0 ]; then
   		stokefiles=($f'A' ${stokefiles[@]})
   	fi
   done
   
   #files=()
   if [ ! ${#stokefiles[@]} -eq 0 ]; then
   	echo 'Make stokes'
   	message=$(sbatch -J 'make_stokes' --mem=5G -t 00:30:00  ~/capo/ael/mk_stokes.py ${stokefiles[@]} --stokes='I')
   	message=($message)
   	id=`echo ${message[3]}`
   	echo $id
   	#Wait for the job to finish
   	while [ `myq | grep $id | wc -l` -ge 1 ]; do
   	    sleep 10
   #	    files=( ${files[@]} $(find . -name "*P" -type 'd' -mmin -0.166)  )   #Find MIRIAD files made less than 10 seconds ago ending in P, append
   	done
   	echo 'mk_stokes completed.'
   fi
   
   for i in "${!files[@]}";
   do
   	files[$i]=${files[$i]}"AP"   #Ensure updated file names are kept.
   done
   
   
   #If files contain the strings "even" or "odd", then split them accordingly into sub-directories.
   #If not, then make separate copies of the files, appending "even" and "odd".
   
   
   #files=$(find . -name "*P" -type 'd' -mtime -10)   #Find MIRIAD files made less than 10 seconds ago ending in P
   
   mkdir -p even
   mkdir -p odd
   cp -u $cal'.py' even/
   cp -u $cal'.py' odd/
   
   fl=1
   for f in ${files[@]}
   do
   	if [[ $f == *"even"* ]]; then
   		cp $f even -ru
   	elif [[ $f == *"odd"* ]]; then
   		cp $f odd -ru
   	else
   	     if [[ $fl == 1 ]]; then
   		     echo 'Warning: Even/odd not found. Making duplicates: ' $
   		     fl=0
   	     fi
                cp $f even/$f -ru
                cp $f odd/$f -ru
   	fi
   done
   
   #!! Skip lst binning. No need for so little data.
   
   #message=$(sbatch -p jpober-test ~/extra_scripts/batch_lstbin_v02.sh $dir)
   #message=$(sbatch ~/extra_scripts/batch_lstbin_v02.sh $calfile)
   
   #~/extra_scripts/batch_lstbin_v02.sh $calfile
   
   #message=($message)
   #id=`echo ${message[3]}`
   #echo $id
   ##Wait for the job to finish
   #while [ `myq | grep $id | wc -l` -ge 1 ]; do
   #    sleep 10
   #done
   
   
   declare -a seps=('0,1' '-1,1' '1,1')
   declare -a evenodd=('even' 'odd')
   for split in "${evenodd[@]}"; do
       cd ${split}
       for sep in "${seps[@]}"; do
           echo $split $sep
           mkdir -p ./'sep'$sep
           bls=$(python ~/capo/pspec_pipeline/getbls.py --sep=${sep} -C $cal ./*AP)
           cd ./'sep'$sep
           python ~/capo/ael/pull_antpols.py -q -a $bls ../*AP
   #       cp -rf -u *.APA ./sep${sep}
   #       rm -rf *.APA
           cd ..
           done
       cd ..
   done
else
   echo 'Skipping to the plotting part'
fi

${SCRIPTSDIR}/mk_psa64_pspec.sh # pspec_psa128.cfg


pwd 
rm $cal.py*
rm even/$cal.py*
rm odd/$cal.py*

cd $olddir
