#! /bin/bash
#SBATCH -t 0:10:00
#SBATCH -n 2                                                                                                                  
#SBATCH --array=0-36:1
###SBATCH --ntasks=1                                                                                                           
#SBATCH --mem=20G                                                                                                              
#SBATCH -J PGInHorizon
####SBATCH -p jpober-test 
#SBATCH --output=/users/jkerriga/brownscripts/jrk/SlurmOut/PGInHorizon_%A_%a.out 
###SBATCH --output=/users/jkerriga/data/jkerriga/PFHDOutput/fhd_%a/FGSub_%A_%a.out
###SBATCH --error=/users/jkerriga/data/jkerriga/PFHDOutput/fhd_%a/FGSub_%A_%a.err 
source activate PAPER
num=$(($SLURM_ARRAY_TASK_ID + 0))
vsname=''
outdir=/users/jkerriga/data/jkerriga/PGProcessed/even/sep0,1

declare -a dirt_list
declare -a model_list

cd $outdir

dirt_list=(lst*uvSBAL)
residual_list=(lst*uvHBAL)

echo "${residual_list[$num]}"

#for num in {0..36};do
echo "${dirt_list[$num]}"
python /users/jkerriga/brownscripts/jrk/ratio.py --dirty "${dirt_list[$num]}" --residual "${residual_list[$num]}"
#done



