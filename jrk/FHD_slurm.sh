#!/bin/bash

####################################################
#
# PIPE_SLURM.SH
#
# Top level script to run a list of observation IDs through FHD (deconvolution or firstpass),
# check the status of resulting FHD outputs, rerun specific observation IDs as necessary
# with new resource allocations, integrate cubes, and generate power spectra through 
# eppsilon.
#
# Required input arguments are obs_file_name (-f /path/to/obsfile) and version
# (-v yourinitials_jackknife_test)
#
# Optional input arguments are: starting_obs (-s 1061311664) which is defaulted to the beginning
# obsid of the specified file, ending_obs (-e 1061323008) which is defaulted to the ending obsid
# of the specified file, outdir (-o /path/to/output/directory) which is defaulted to 
# /nfs/mwa-09/r1/djc/EoR2013/Aug23, wallclock_time (-w 08:00:00) which is defaulted to 
# 4 hours for a typical firstpass run, ncores (-n 10) which is defaulted to 10 for a typical IDL
# job, mem (-m 4G) which is defaulted to 4 Gigabytes per slot for a typical firstpass run, and
# thresh (-t 1) which is defaulted to 1 to tell the code to not look for a threshold from wedge
# statistics.
#
# WARNING!
# Terminal will hang as it waits for jobs to finish, and closing the termianal will kill any 
# remaining jobs! To run in the background, run: 
# nohup ./pipe_slurm.sh -f /path/to/obsfile -v yourinitials_jackknife_test > /path/to/your/output/log/file.txt &
#
####################################################

#Clear input parameters
unset obs_file_name
unset starting_obs
unset ending_obs
unset outdir
unset version
unset resubmit_list
unset resubmit_index

#######Gathering the input arguments and applying defaults if necessary

#Parse flags for inputs
#while getopts ":f:s:e:o:v:p:w:n:m:t:" option
while getopts ":f:s:e:o:v:w:n:m:t:" option
do
   case $option in
	f) obs_file_name="$OPTARG";;	#text file of observation id's
	s) starting_obs=$OPTARG;;	#starting observation in text file for choosing a range
	e) ending_obs=$OPTARG;;		#ending observation in text file for choosing a range
        o) outdir=$OPTARG;;		#output directory for FHD output folder
        v) version=$OPTARG;;		#FHD folder name and case for eor_firstpass_versions
					#Example: nb_foo creates folder named fhd_nb_foo
	w) wallclock_time=$OPTARG;;	#Time for execution in slurm
	n) ncores=$OPTARG;;		#Number of cores for slurm
	m) mem=$OPTARG;;		#Memory per node for slurm
	t) thresh=$OPTARG;;		#Wedge threshold to use to determine whether or not to run
        p) partition=$OPTARG;;
       \?) echo "Unknown option: Accepted flags are -f (obs_file_name), -s (starting_obs), -e (ending obs), -o (output directory), "
	    echo "-v (version input for FHD), -w (wallclock timein slurm), -n (number of cores to use),"
	    echo "and -m (memory per core for slurm)." 
	    exit 1;;
	:) echo "Missing option argument for input flag"
	   exit 1;;
   esac
done

#Manual shift to the next flag.
shift $(($OPTIND - 1))

#Specify the FHD file path that is used in IDL (generally specified in idl_startup)
FHDpath=$(idl -e 'print,rootdir("fhd")') ### NOTE this only works if idlstartup doesn't have any print statements (e.g. healpix check)
if [ -z ${wallclock_time} ]; then
    wallclock_time=5:00:00
fi
#Set typical nodes needed for standard FHD firstpass if not set.                                                                                                     
if [ -z ${ncores} ]; then
    ncores=20
fi
#Set typical memory needed for standard FHD firstpass if not set.                                                                                                    
if [ -z ${mem} ]; then
    mem=80G
fi
if [ -z ${partition} ]; then
    partition='default-batch'
fi
if [ -z ${thresh} ]; then
    # if thresh is not set, set it to -1 which will cause it to not check for a window power                                                                         
    thresh=-1
fi

if [ -z ${outdir} ]
then
    outdir=/users/jkerriga/data/jkerriga/2DayOutput
    echo Using default output directory: $outdir
else
    #strip the last / if present in output directory filepath
    outdir=${outdir%}
    echo Using output directory: $outdir
fi

mkdir -p ${outdir}/fhd_${version}
mkdir -p ${outdir}/fhd_${version}/grid_out
echo Output located at ${outdir}/fhd_${version}

#Read the obs file and put into an array, skipping blank lines if they exist
#i=0
#while read line
#do
#   if [ ! -z "$line" ]; then
#      obs_id_array[$i]=$line
#      i=$((i + 1))
#   fi
#done < "$obs_file_name"
#obs_list=$obs_file_name
#obs_list='Pzen.2456242.45916.uvcRREcACOTUc'
#obs_list='Pzen.2456242.30605.uvcRREcACOTUc.uvfits'
obs_list='Pzen.2456242.48003.uvcRREcACOTUc.uvfits'
#obs_list='PPzen.2456242.41044.uvcRREcACOTUcSB.uvfits'
#obs_list='zen.2456242.29909.uvcRREcACOM.uvfits'
#obs_list='Szen.2456327.34089.uvcRREcACOM.uvfits'
export $version
message=$(sbatch --mem=$mem -t ${wallclock_time} -n ${ncores} -p ${partition} --export=ncores=$ncores,outdir=$outdir/fhd_${version},version=$version,thresh=$thresh -o ${outdir}/fhd_${version}/grid_out/firstpass-%A_%a.out -e ${outdir}/fhd_${version}/grid_out/firstpass-%A_%a.err ~/brownscripts/jrk/eor_FHD_slurm_job.sh ${obs_list})

#echo $message

#Run the command
message=($message)

echo ${message[@]}

#Gather the job id from the job for later use

##
id=`echo ${message[3]}`
echo $id



########End of submitting the firstpass job and waiting for output


while [ `myq | grep $id | wc -l` -ge 1 ]; do
    sleep 10
done


#Check to see if Healpix cubes exist for all obsids
i=0
rerun_flag=0
for obs_id in "${obs_id_array[@]}"; do
    i=$((i + 1))
    # Check to see if 4 files (even/odd, XX/YY) return from listing for that obsid
    if ! ls -1 ${outdir}/fhd/Healpix/${obs_id}*cube* 2>/dev/null | wc -l | grep 4 -q; then
	echo Observation $obs_id is missing one or more Healpix cubes
        rerun_flag=1
        [[ $resubmit_list =~ $x ]] || resubmit_list+=($obs_id)
        [[ $resubmit_index =~ $i ]] || resubmit_index+=($i)
    fi

done


