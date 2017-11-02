#! /bin/bash  

declare -a seps=('0,1' '-1,1' '1,1')
baselines=()
#cd /users/jkerriga/data/jkerriga/PSA64FHDLST
for sep in "${seps[@]}"; do
    bls=$(python ~/capo/pspec_pipeline/getbls.py --sep=${sep} -C psa6240_FHD ~/data/jkerriga/PSA64FHD/zen.2456242.30605.uvcRREcACOTUc)
    baselines=($bls','$baselines)
done
echo $baselines
python ~/capo/dcj/scripts/pull_antpols.py -a $baselines $* #Pzen.2456*.*.*c[HS]P

#mv Pzen.2456*[02468].*.uvcRREcACOTUc[SH]PA ./even/


    
     