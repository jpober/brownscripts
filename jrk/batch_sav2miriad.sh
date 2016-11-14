#! /bin/bash
#SBATCH -t 0:25:00
#SBATCH -n 2                                                                                                                  
#SBATCH --array=0-999:1
###SBATCH --ntasks=1                                                                                                           
#SBATCH --mem=20G                                                                                                              
#SBATCH -J Sav2Miriad
####SBATCH -p jpober-test 
#SBATCH --output=/users/jkerriga/brownscripts/jrk/SlurmOut/Sav2Mir_%A_%a.out 
###SBATCH --output=/users/jkerriga/data/jkerriga/PFHDOutput/fhd_%a/FGSub_%A_%a.out
###SBATCH --error=/users/jkerriga/data/jkerriga/PFHDOutput/fhd_%a/FGSub_%A_%a.err 
source activate PAPER
#2300
version=$(($SLURM_ARRAY_TASK_ID + 999))
vsname=''
outdir=/users/jkerriga/data/jkerriga/AnalysisOutput

obs_list=($(cat Analysis_obsfits.txt))

echo ${obs_list[$version]}
#obs_id=${obs_list[$version]}

cd ${outdir}/fhd_${vsname}${version}/vis_data
python ~/brownscripts/jrk/sav2miriad.py ./*sav ../metadata/*
#python ~/brownscripts/jrk/mk_stokes.py --stokes='I' Pzen*[SH]
#python ~/brownscripts/jrk/WideBandFilter.py Pzen*[SH]P
#~/brownscripts/jrk/reduce_seps.sh

#python ~/capo/pspec_pipeline/pspec_prep.py -C psa6240_FHD --model -a cross --nogain --nophs --clean=1e-7 --horizon=0 --window='none' Pzen*uvcRREcACOTUcMP
python ~/brownscripts/jrk/WideBandFilter.py Pzen*uvcRREcACOTUcMP

#rm -r *SP
#uv_addsub.py --sub *HP *MPF
~/brownscripts/jrk/reduce_seps.sh *HP *MPF *MP
python ~/brownscripts/jrk/uv_sub.py *HPA *MPFA -s SPA
#~/brownscripts/jrk/batch_submodel.sh
#cp -rf *[HS]P ../../../PaperAnalysis/AllSeps/
if [ $(du -sh *SPA|wc -l) = 1 ];then
    cp -rf Pzen*HPA Pzen*SPA ../../../Analysis-10LST/
    cp -r Pzen*MPA Pzen*MPFA ../../../Analysis-10Models/
fi
rm -r *Uc[HSM]*


