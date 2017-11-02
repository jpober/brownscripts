#! /bin/bash
#SBATCH -t 2:00:00
#SBATCH -n 4
#SBATCH --mem=40G
####SBATCH -p jpober-test
#####SBATCH --array=0-95:1
#SBATCH --array=0-95:1
#SBATCH -J LSTBin
source activate PAPER
#day=('even' 'odd')
LSTS=`python -c "import numpy as n; hr = 1.0; print ' '.join(['%f_%f' % (d,d+hr/4.) for d in n.arange(0,24*hr,hr/4.)])"`
echo $LSTS

MY_LSTS=($(~/capo/dcj/scripts/pull_args.py $LSTS))
CALFILE=hsa7458_v000

EXT=zen*uvc

#cp psa6240_FHD.py /users/jkerriga/data/jkerriga/PSA64FHDLST/${day}/
#for i in ${day[@]};do
cd /users/jkerriga/data/shared/hera19season 
#cd /users/jkerriga/data/jkerriga/HERA2
python ~/capo/zsa/scripts/lstbin_v02.py -a cross -C ${CALFILE} -s Sun --lst_res=10.737 --lst_rng="${MY_LSTS[$SLURM_ARRAY_TASK_ID]}" --tfile=601.272 --altmax=0 --stats='cnt' ${EXT} 
