#! /bin/bash  

declare -a seps=('0,1' '-1,1' '1,1')
declare -a evenodd=('even' 'odd')
for split in "${evenodd[@]}"; do
    cd ~/data/jkerriga/UnphasedConFit/${split}
    for sep in "${seps[@]}"; do
	echo $split $sep
	mkdir ./'sep'$sep
	bls=$(python ~/capo/pspec_pipeline/getbls.py --sep=${sep} -C psa6240_FHD ~/data/jkerriga/UnphasedConFit/lst.2456242.30605.30605.uvS)
	python ~/capo/scripts/pull_antpols.py -a $bls ~/data/jkerriga/UnphasedConFit/${split}/lst*uvS
	cp -rf *A ./sep${sep}
	rm -r *A
	done
done
    
     