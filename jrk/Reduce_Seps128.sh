#! /bin/bash  

declare -a seps=('0,1' '-1,1' '1,1')
baselines=()
#cd /users/jkerriga/data/jkerriga/PSA64FHDLST
for sep in "${seps[@]}"; do
    bls=$(python ~/capo/pspec_pipeline/getbls.py --sep=${sep} -C psa6240_FHD ~/data/shared/PAPER128_Simulations/PAPER128_0)
    baselines=($bls','$baselines)
done
echo $baselines
python ~/capo/dcj/scripts/pull_antpols.py -a $baselines PAPER128_*

#mv Pzen.2456*[02468].*.uvcRREcACOTUc[SH]PA ./even/


    
     