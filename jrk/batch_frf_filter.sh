#!/bin/bash

seps='sep0,1 sep1,1 sep-1,1'
cal=psa6240_FHD
out=/users/jkerriga/data/jkerriga/PaperData
appelation=uvHA   #uvcRREcACOTUcHPA
chan='101' #was 101

scriptsdir=~/scratch/capo_jrk/mjk/scripts

bl_scale='1.2' #was 1.2
fr_width='1.3' #was 1.3
maxfr='0.0013' #was 0.0013
alietal=true

declare -A ants
ants[sep0,1]=41_49 #47_49
ants[sep1,1]=1_3 #1_3
ants[sep-1,1]=1_48 #1_48

days='evenG oddG'

printf 'This is Batch Fringe Rate Filter using:\n'
printf 'calfile: %s \n' $cal
printf 'bl_scale: %s\n' $bl_scale
printf 'fr_width_scale: %s\n' $fr_width
printf 'maxfr: %s\n' $maxfr
sleep 1.5

if [ ! -d $out   ]; then
    printf 'This output director appears to be new\n'
    printf 'Creating Directoires now\n'
    for path in $paths; do
        for sep in $seps; do
            mkdir -p ${out}/${path}/${sep}
        done
   done
fi

for path in $days; do
    for sep in $seps; do
        printf 'Checking for Old FRF Data\n'
        outpath=$out/$path
        if [[ $(ls -d $outpath/$sep/*.${appelation}L) ]] ; then
            printf 'Found %s folders\n' $(ls -d $outpath/$sep/*.${appelation}L| wc -l)
            printf 'First deleting old FRF Filtered Data \n'
            sleep 1
            printf 'I wil delete: \n'
            sleep 1
            printf '%s \n' $(ls -d $outpath/$sep/*.${appelation}L)
            read -p "Are you sure you want to delete these? [y/N] " -n 1 -r
            printf '\n'
            if [[ $REPLY =~ ^[Yy]$ ]]
                then
                    rm -rf $(ls -d $outpath/$sep/*.${appelation}L)
                    printf 'Files are deleted\nContinuting to Fringe Rate Filter\n'
                else
                    printf 'Nothing is deleted\ncontinuing\n'
                    continue
            fi
        fi

        files=$(ls -d $outpath/$sep/*.${appelation})
        printf 'Filtering %s by selecting ant %s \n' $sep ${ants[$sep]}

        if ${alietal}; then

            printf '%s/frf_filter.py -C %s  --alietal -a %s -c %s --bl_scale %s --fr_width_scale %s --outpath=%s/' \
            $scriptsdir  $cal ${ants[$sep]} $chan $bl_scale $fr_width ${out} 
            printf '%s\n' $path/$sep
            "${scriptsdir}/frf_filter.py" -C $cal -a ${ants[$sep]} -C $cal --alietal --bl_scale $bl_scale \
                --fr_width_scale $fr_width $files --outpath='' -c $chan --maxfr=${maxfr} --pol=I 

        else
            printf '%s/frf_filter.py -C %s   -a %s -c %s --bl_scale=%s --fr_width_scale=%s --outpath=%s/' $scriptsdir $cal ${ants[$sep]} $chan $bl_scale $fr_width ${out} 
            printf '%s\n' $path/$sep
            "${scriptsdir}/frf_filter.py" -C $cal -a ${ants[$sep]} -C $cal  $files --outpath='' -c $chan --pol=I  --maxfr=${maxfr}  --bl_scale $bl_scale --fr_width_scale=${fr_width}

        fi
    done
done
