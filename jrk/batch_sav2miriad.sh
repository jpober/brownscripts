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
version=$(($SLURM_ARRAY_TASK_ID + 0))
vsname=''
outdir=/users/jkerriga/data/jkerriga/AnalysisOutput

#obs_list=($(cat Analysis_obsfits.txt))

#echo ${obs_list[$version]}
#obs_id=${obs_list[$version]}

cd ${outdir}/fhd_${vsname}${version}/vis_data
python ~/brownscripts/jrk/sav2miriad.py ./*sav ../metadata/*
python ~/brownscripts/jrk/ModelDelayFilter.py Pzen*uvcRREcACOTUcMP

#rm -r *SP
#uv_addsub.py --sub *HP *MPF

~/brownscripts/jrk/reduce_seps.sh *HP *MPF *MP
python ~/brownscripts/jrk/uv_sub.py --dirty *HPA --model *MPFA -s SPA

#~/brownscripts/jrk/batch_submodel.sh
#cp -rf *[HS]P ../../../PaperAnalysis/AllSeps/
if [ $(du -sh *SP|wc -l) = 1 ];then
    cp -rf Pzen*HP Pzen*SP ../../../AllSepsLST/
    #cp -r Pzen*MPA Pzen*MPFA ../../../PGModels/
fi
#rm -r *Uc[HSM]*


