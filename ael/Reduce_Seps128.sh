#! /bin/bash  

#Arguments in = List of Miriad files

declare -a seps=('0,1' '-1,1' '1,1')
baselines=()
#cd /users/jkerriga/data/jkerriga/PSA64FHDLST
for sep in "${seps[@]}"; do
    bls=$(python ~/capo/pspec_pipeline/getbls.py --sep=${sep} -C psa6622_v003_paper128 $1)
    baselines=($bls','$baselines)
done
echo $baselines
#Use 
python ~/capo/dcj/scripts/pull_antpols.py -a $baselines $@

#mv Pzen.2456*[02468].*.uvcRREcACOTUc[SH]PA ./even/
