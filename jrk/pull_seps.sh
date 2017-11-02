#! /bin/bash  

declare -a seps=('0,1' '-1,1' '1,1')
declare -a evenodd=('evenIH' 'oddIH')
dir=PaperData
suffix=uvHBMx
for split in "${evenodd[@]}"; do
    cd ~/data/jkerriga/${dir}/${split}
    for sep in "${seps[@]}"; do
	echo $split $sep
	mkdir ./'sep'$sep
	bls=$(python /users/jkerriga/scratch/capo_jrk/pspec_pipeline/getbls.py --sep=${sep} -C psa6240_FHD ~/data/jkerriga/PSA64LST/zen.2456245.57048.uvcRREcACOTUcP)
	python ~/capo/dcj/scripts/pull_antpols.py -a $bls ~/data/jkerriga/${dir}/${split}/lst*${suffix}
	cp -rf lst*${suffix}A ./sep${sep}
	rm -rf lst*${suffix}A
	done
done
    
     