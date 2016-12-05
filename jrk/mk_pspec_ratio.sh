#! /bin/bash
export PYTHONPATH='.':$PYTHONPATH
#PREFIX="OneDayFG"
#
##chans=`python -c "print ' '.join(['%d_%d'%(i,i+39) for i in range(10,150,1)])"`
#pols='I Q U V'
#seps='0_16 1_16 0_17'
#chans='110_149'
#RA="1:01_9:00"
#NBOOT=20
#
##DATAPATH=fringe_hor_v006
#SCRIPTSDIR=~/src/capo/pspec_pipeline
#cal="psa898_v003"
#PWD=`pwd`
#DATAPATH="${PWD}/typical_day/*FRXS"

(
echo using config $*
. $*

#set defaults to parameters which might not be set
if [[ ! -n ${WINDOW} ]]; then export WINDOW="none"; fi
#for posterity Print the cal file we are using

cov_array=($covs)
chan_array=($chans)
if [ -z "$covs" ]
then
    declare -A scales
    for (( i=0; i<${#chan_array[@]}; ++i )); do
        scales[${chan_array[$i]}]=1
    done
else
    declare -A scales
    for (( i=0; i<${#chan_array[@]}; ++i )); do
        scales[${chan_array[$i]}]=${cov_array[$i]}
    done
fi


pywhich $cal


threadcount=`python -c "c=map(len,['${pols}'.split(),'${chans}'.split(),'${seps}'.split()]);print c[0]*c[1]*c[2]"`
echo Running $threadcount pspecs
PIDS=""

#FILES=`lst_select.py -C ${cal} --ra=${RA} ${DATAPATH}`
test -e ${PREFIX} || mkdir $PREFIX
for chan in $chans; do
    chandir=${PREFIX}/${chan}
    test -e ${chandir} || mkdir ${chandir}
    for pol in $pols; do
        echo "Starting work on ${pol}" 
        poldir=${chandir}/${pol}
        test -e ${poldir} || mkdir ${poldir}
        if [ ! -e ${poldir}/pspec_${PREFIX}_${chan}_${pol}.png ]; then
            for sep in $seps; do
                sepdir=${poldir}/${sep}
		oddLSTS=(${ODD_DATAPATH}${sep}/lst.*.*.*.${SUFFIX})
		evenLSTS=(${EVEN_DATAPATH}${sep}/lst.*.*.*.${SUFFIX})
		oddLSTS2=(${ODD_DATAPATH}${sep}/lst.*.*.*.${SUFFIX2})
		evenLSTS2=(${EVEN_DATAPATH}${sep}/lst.*.*.*.${SUFFIX2})
                EVEN_FILES="${evenLSTS[@]:8:28}" #${EVEN_DATAPATH}${sep}/ #lst.*.[345]*.*.${SUFFIX} #uvHBFAL was 3:15
                ODD_FILES="${oddLSTS[@]:8:28}" #${ODD_DATAPATH}${sep}/ #lst.*.[345]*.*.${SUFFIX} #uvHBFAL
		EVEN_FILES2="${evenLSTS2[@]:8:28}"
		ODD_FILES2="${oddLSTS2[@]:8:28}"
		echo $EVEN_FILES
                test -e ${sepdir} || mkdir ${sepdir}
                LOGFILE=`pwd`/${PREFIX}/${chan}_${pol}_${sep}.log
                echo this is mk_psa64_pspec.sh with  |tee  ${LOGFILE}
                echo experiment: ${PREFIX}|tee -a ${LOGFILE}
                echo channels: ${chan}|tee -a ${LOGFILE}
                echo polarization: ${pol}|tee -a ${LOGFILE}
                echo separation: ${sep}|tee -a ${LOGFILE}
                echo `date` | tee -a ${LOGFILE}

                #ANTS=`grid2ant.py -C ${cal} --seps="${sep}"`
                ANTS='cross'
                echo python ${SCRIPTSDIR}/pspec_cov_v002.py -C ${cal} \
                     -b ${NBOOT} -a ${ANTS} -c ${chan} -p ${pol}\
                      --window=${WINDOW}  ${NOPROJ} --output=${sepdir} \
                       ${EVEN_FILES} ${ODD_FILES} 
                
                python ${SCRIPTSDIR}/jrk_pspec_cov_v003.py -C ${cal} -b ${NBOOT} \
                    -a ${ANTS} -c ${chan} -p ${pol} --window=${WINDOW} \
                      ${NOPROJ} --output=${sepdir} \
                      --de "${EVEN_FILES}" --do "${ODD_FILES}" --re "${EVEN_FILES2}" --ro "${ODD_FILES2}" \
                     | tee -a ${LOGFILE}

                
                echo beginning bootstrap: `date` | tee -a ${LOGFILE} 
                ${SCRIPTSDIR}/jrk_pspec_cov_boot_v002.py --identity ${sepdir}/pspec_boot*npz | tee -a ${LOGFILE} 
                echo complete! `date`| tee -a ${LOGFILE} 
                mv pspec.npz ${sepdir}/
                PIDS="${PIDS} "$!
            done
        fi
    done
done

echo waiting on `python -c "print len('${PIDS}'.split())"` power spectra threads ${PIDS} 
wait $PIDS
echo power spectrum complete






echo averaging power spectra for pols/channels
for chan in $chans; do
    chandir=${PREFIX}/${chan}
    for pol in $pols; do
        echo "Generating plots for ${chan}: ${pol}"
        echo "Multiplying pspec by factor ${scales[${chan}]} for Cov"
        poldir=${chandir}/${pol}
        #PLOT
        ${SCRIPTSDIR}/plot_pk_k3pk_zsa_2.py ${poldir}/*/pspec.npz --cov=${scales[$chan]} #--flux
        mv pspec_pk_k3pk.npz pspec_${PREFIX}_${chan}_${pol}.npz
        mv pspec.png pspec_${PREFIX}_${chan}_${pol}.png
        mv posterior.png posterior_${PREFIX}_${chan}_${pol}.png
        mv posterior.txt posterior_${PREFIX}_${chan}_${pol}.txt
        cp  pspec_${PREFIX}_${chan}_${pol}.png ${poldir}/
        cp  posterior_${PREFIX}_${chan}_${pol}.png ${poldir}/
        cp  pspec_${PREFIX}_${chan}_${pol}.npz ${poldir}/
        cp  posterior_${PREFIX}_${chan}_${pol}.txt ${poldir}/
    done
done
)
