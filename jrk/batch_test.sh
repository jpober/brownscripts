#! /bin/bash
#SBATCH -t 5:00:00
#SBATCH -n 8
#SBATCH --mem=50G
#SBATCH -p jpober-test

LSTS=(`python -c "import numpy as n; hr = 1.0; print ' '.join(['%f_%f' % (d,d+hr/4.) for d in n.arange(0,23.999*hr,hr/4.)])"`)

#LSTS=`python -c "import numpy as n; hr = 1.0; print ' '.join(['%f_%f' % (d,d+hr/4.) for d in n.arange(6,9*hr,hr/4.)])"`
#LSTS=`python -c "import numpy as n; hr = 1.0; print ' '.join(['%f_%f' % (d,d+hr/4.) for d in n.arange(5.25,9.5*hr,hr/4.)])"`
#LSTS=9.250000_9.500000
MY_LSTS=($(python ~/capo/dcj/scripts/pull_args.py $LSTS))
CALFILE=psa6240_FHD
EXT=uvcRREcACOTUcP
trap "exit" INT
days=('even' 'odd')
echo ${MY_LSTS[0]}
echo ${days[1]}