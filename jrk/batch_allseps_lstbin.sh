#! /bin/bash
#SBATCH -t 2:00:00
#SBATCH -n 2
#SBATCH --mem=20G
#####SBATCH -p jpober-test
#####SBATCH --array=0-95:1
#SBATCH --array=0-36:1
#SBATCH -J LSTBin

source activate PAPER
LSTS=`python -c "import numpy as n; hr = 1.0; print ' '.join(['%f_%f' % (d,d+hr/4.) for d in n.arange(0,24*hr,hr/4.)])"`

#MY_LSTS=`~/capo/dcj/scripts/pull_args.py $LSTS`
MY_LSTS=($(~/capo/dcj/scripts/pull_args.py $LSTS))
CALFILE=psa6240_FHD
#EXT=*uvcRREcACOTUcSPA
EXT=*HP
#obs=$(cat ~/data/jkerriga/FGSubLST/${day}Obs.txt)
#obsfixed=()
#for j in ${obs[*]};do
#    echo ${j%txt}$EXT
#    obsfixed+=(${j%txt}$EXT)
#done
evenodd='even odd'
#echo "${obsfixed[@]}"
#EXT=*uvcRREcACOTUcSPAB #remember to change file suffix too
#cp psa6240_FHD.py /users/jkerriga/data/jkerriga/PSA64FHDLST/${day}/
for day in ${evenodd}; do
    cd /users/jkerriga/data/jkerriga/AllSepsLST/${day} #/${i}
    echo ${day}
    echo "${MY_LSTS[$SLURM_ARRAY_TASK_ID]}"
    python ~/capo/scripts/lstbin_v02.py -a cross -C ${CALFILE} -s Sun --lst_res=42.95 --lst_rng="${MY_LSTS[$SLURM_ARRAY_TASK_ID]}" --tfile=600 --altmax=0 --stats='cnt' --median --nsig=3 ${EXT} #"${obsGood[@]}" #"${obsBad[@]}" #${EXT}
done