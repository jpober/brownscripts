#! /bin/bash                                                                                                                                     
#SBATCH -t 4:00:00                                                                                                                               
#SBATCH --array=242-249:1
#SBATCH -n 10                                                                                                                                   
#SBATCH --mem=10G                                                                                                                                
####SBATCH -p jpober-test                                                                                                                           
####SBATCH --output=/users/jkerriga/brownscripts/jrk/SlurmOut/FGSub_%A_%a.out  
source activate PAPER

version=$(($SLURM_ARRAY_TASK_ID + 0))
outdir=/users/jkerriga/data/jkerriga/NewSubLST

module load ghostscript
module load imagemagick/6.6.4
module load git/2.2.1
cd ${outdir}
#python ~/capo/zsa/scripts/pspec_prep.py -C psa6240_FHD -a cross --nogain --nophs --clean=1e-4 --horizon=15 --window='blackman-harris'  Pzen.*${version}.*.*[SH]PA

python ~/brownscripts/jrk/WideBandFilter.py Pzen.*${version}.*.*HPA*
#python ~/capo/dcj/scripts/xrfi_simple.py -n 3 Pzen.*${version}.*.*[SH]PAF

#python ~/capo/pspec_pipeline/pspec_prep.py -C psa6240_FHD -a cross --window=blackman-harris --nogain --nophs --clean=1e-9 --horizon=15 Pzen.*.${version}*.*[SH]P
