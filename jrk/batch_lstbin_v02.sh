#! /bin/bash
#SBATCH -t 4:00:00
#SBATCH -n 2
#SBATCH --mem=20G
#####SBATCH --array=0-95:1
#SBATCH --array=0-36:1
#SBATCH -J LSTBin
#SBATCH --qos=jpober-condo

source activate PAPER

## even/odd day sep names ##
day=('evenIH' 'oddIH')
#day=('Bin2' 'Bin8' 'Bin16' 'Bin34')
## get list of lsts ##
LSTS=`python -c "import numpy as n; hr = 1.0; print ' '.join(['%f_%f' % (d,d+hr/4.) for d in n.arange(0,24*hr,hr/4.)])"`
MY_LSTS=($(~/capo/dcj/scripts/pull_args.py $LSTS))

## cal file used for obs ##
CALFILE=psa6240_FHD

## extension of obs being LST binned ##
EXT=*SPABR

## directory of even/odd split folders ##
obsDir=/users/jkerriga/data/jkerriga/PaperData

## loops over even/odd days ##
for i in ${day[@]};do
    cd ${obsDir}/${i}
    echo $i
    python ~/capo/scripts/lstbin_v02.py -a cross -C ${CALFILE} -s Sun --lst_res=42.95 --lst_rng="${MY_LSTS[$SLURM_ARRAY_TASK_ID]}" --tfile=600 --altmax=0 --stats='cnt' --median --nsig=3 ${EXT}

done

