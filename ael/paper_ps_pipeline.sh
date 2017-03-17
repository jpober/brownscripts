#!/bin/bash

## A wrapper to run the PAPER power spectrum pipeline on FHD simulated data.

#Parse input arguments. All other arguments passed are the MIRIAD files to process.
touch tmp
skiptoplot=0
reset=0


while [[ $@ ]]
do
	OPTIND=1
	if [[ $1 =~ ^- ]]
	then
	    getopts ":c:sr" option
	    case $option in
		c) configfile="$OPTARG"
		   shift ;;
		s) skiptoplot=1 ;;
		r) reset=1
	    esac
	else
	    echo $1 >> tmp
	fi
	shift
done


fileItemString=$(cat tmp |tr "\n" " ")
infiles=( $fileItemString )
rm tmp

if [ $reset -eq 1 ]; then
   for f in ${infiles[@]};
   do
   	mv $f'.' $f   #Optionally, remove the dots from original filenames
   done
   exit
fi
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
  

if [ -z ${auto_redundancy} ]; then
	auto_redundancy=0
fi

if [ ! -z ${non_redundant} ]; then
	## The array is non-redundant, so run mk_nonred_pspec.sh instead.
	skiptoplot=1
	
fi

if [ -z ${length} ]; then
	length=0
fi


if [ $auto_redundancy == 1 ]; then
	s=$(python ~/capo/ael/build_redundant.py --count -C $cal --length=$length --min=8)	#Return the number of redundant groups
	s=$(( s - 1 ))
	seps=( $(seq 0 $s) )
	echo "Seps: " ${seps[@]}
#	seps=( "${seps[@]/#/sep}" )
	declare -a seps=${seps[@]}
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
   		mv $f $f'.' 2>/dev/null  #Important --> Needs a suffix
   		files=( ${files[@]} $f'.' )
   	else
   		files=( ${files[@]} $f )
   	fi
   done
   echo 'Reducing files'
   reducfiles=()
   echo ${files[@]}
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

	if [ $auto_redundancy == 1 ]; then
	# Build seps list from redundancies
		echo "Auto Redundancies"
		baselines=$(python ~/capo/ael/build_redundant.py --min=8 --length=$length --flatten -C $cal --save=$cal"_reds")    #Return a flattened list of all redundant baselines
		#s=$(python ~/capo/ael/build_redundant.py --count -C $cal --restore=$cal_"reds.npz")	#Return the number of redundant groups
		#declare -a seps=$(seq 0 $s)
	else
		baselines=()
		#cd /users/jkerriga/data/jkerriga/PSA64FHDLST
		for sep in "${seps[@]}"; do
		    bls=$(python ~/capo/pspec_pipeline/getbls.py --sep="${sep}" -C $cal ${reducfiles[1]})
		    baselines=($bls','$baselines)
	done
	fi

	echo $baselines
	#Use 
	python ~/capo/dcj/scripts/pull_antpols.py -a $baselines ${reducfiles[@]}


#   	Reduce_Seps128.sh ${reducfiles[@]}
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
   	#message=$(sbatch -J 'make_stokes' --mem=5G -t 00:30:00  ~/capo/ael/mk_stokes.py ${stokefiles[@]} --stokes='I')
   	~/capo/ael/mk_stokes.py ${stokefiles[@]} --stokes='I'
   	#message=($message)
   	#id=`echo ${message[3]}`
   	#echo $id
   	#Wait for the job to finish
   	#while [ `myq | grep $id | wc -l` -ge 1 ]; do
   	    sleep 10
   #	    files=( ${files[@]} $(find . -name "*P" -type 'd' -mmin -0.166)  )   #Find MIRIAD files made less than 10 seconds ago ending in P, append
   	#done
   	echo 'mk_stokes completed.'
   fi
   
   for i in "${!files[@]}";
   do
   	files[$i]=${files[$i]}"AP"   #Ensure updated file names are kept.
   done
   
   
   #If files contain the strings "even" or "odd", then split them accordingly into sub-directories.
   #If not, then make separate copies of the files, appending "even" and "odd".
   
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
   
   
   declare -a evenodd=('even' 'odd')
   for split in "${evenodd[@]}"; do
       cd ${split}
       for sep in "${seps[@]}"; do
           echo $split $sep
           #mkdir -p ./'sep'$sep
           mkdir -p ./$sep
	   if [ $auto_redundancy == 1 ]; then
		bls=$(python ~/capo/ael/build_redundant.py --length=$length --min=8 --sep=${sep} -C $cal --restore=$cal'_reds.npz')	# Return redundant group
	   else
           	bls=$(python ~/capo/pspec_pipeline/getbls.py --sep=${sep} -C $cal ./*AP)
	   fi
           cd ./$sep
           python ~/capo/ael/pull_antpols.py -q -a $bls ../*AP
   #       cp -rf -u *.APA ./sep${sep}
   #       rm -rf *.APA
           cd ..
           done
       cd ..
   done
else
   echo 'Skipping file setup'
fi


if [ ! -z ${non_redundant} ]; then
   ## The array is non-redundant, so run mk_nonred_pspec.sh instead.
   files=()
   for f in ${infiles[@]};
   do
   #	count=$(python -c "print \""$f"\".count('.')")
   #	if [ $count -eq 0 ]; then
   	if [[ ! $f == *'.A'* ]]; then
   		mv $f $f'.A' 2>/dev/null  #Important --> Needs a suffix
   		files=( ${files[@]} $f'.A' )
   	else
   		files=( ${files[@]} $f )
   	fi
   done

   for f in ${files[@]};
   do
        f=$(basename $f)
        tmp=( $(find . -name $f'P' -type 'd') )
        if [ ${#tmp[@]} -eq 0 ]; then
                stokefiles=($f ${stokefiles[@]})
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
   #        files=( ${files[@]} $(find . -name "*P" -type 'd' -mmin -0.166)  )   #Find MIRIAD files made less than 10 seconds ago ending in P, append
        done
        echo 'mk_stokes completed.'
   fi

   for i in "${!files[@]}";
   do
        files[$i]=${files[$i]}"P"   #Ensure updated file names are kept.
   done


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
                cp $f even/$f'A' -ru
                cp $f odd/$f'A' -ru
   	fi
   done
	. ${SCRIPTSDIR}/mk_nonred_pspec.sh
else
	# File setup already done, so just run the power spectrum
	. ${SCRIPTSDIR}/mk_psa64_pspec.sh # pspec_psa128.cfg
fi



pwd 
#rm $cal.py*
rm even/$cal.py*
rm odd/$cal.py*

cd $olddir
