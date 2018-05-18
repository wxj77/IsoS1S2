#!/bin/bash
#SBATCH -J LUX_Iso_R3
#SBATCH -n 1
#SBATCH -c 24
#SBATCH --mem=48gb 
#SBATCH --mail-type=all    
#SBATCH -t 2-00:00:00
#SBATCH -p standard
#SBATCH -o /scratch/dkhaitan/IsolatedS1S2/salt_output/matlab_real1.log
#SBATCH -a 1-605
module load matlab
matlab -nodesktop -r "Level1_r2_sbatch_2($SLURM_ARRAY_TASK_ID)"  >> /scratch/dkhaitan/IsolatedS1S2/salt_output/log_real1.txt
