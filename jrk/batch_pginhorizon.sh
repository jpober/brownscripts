#! /bin/bash
#SBATCH -t 0:10:00
#SBATCH -n 2                                                                                                                  
#SBATCH --array=0-999:1
###SBATCH --ntasks=1                                                                                                           
#SBATCH --mem=20G                                                                                                              
#SBATCH -J PGInHorizon
####SBATCH -p jpober-test 
#SBATCH --output=/users/jkerriga/brownscripts/jrk/SlurmOut/PGInHorizon_%A_%a.out 
###SBATCH --output=/users/jkerriga/data/jkerriga/PFHDOutput/fhd_%a/FGSub_%A_%a.out
###SBATCH --error=/users/jkerriga/data/jkerriga/PFHDOutput/fhd_%a/FGSub_%A_%a.err 
source activate PAPER
num=$(($SLURM_ARRAY_TASK_ID + 3999))
vsname=''
outdir=/users/jkerriga/data/jkerriga/PGBH
modeldir=/users/jkerriga/data/jkerriga/PGModels
declare -a dirt_list
declare -a model_list

#cd $outdir
#dirt_list=(Pzen*HPA)
cd $modeldir
model_list=(Pzen*MPA)

echo "${model_list[$num]}"
# "${dirt_list[$num]}"

python /users/jkerriga/brownscripts/jrk/ModelDelayFilter.py "${model_list[$num]}"

cp -r "${model_list[$num]}""F" $outdir
rm -r "${model_list[$num]}""F"
cd $outdir
dirt_list=(Pzen*HPA)
python /users/jkerriga/brownscripts/jrk/uv_sub.py --dirty "${dirt_list[$num]}" --model "${model_list[$num]}""F" -s SPA



