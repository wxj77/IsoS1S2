#!/bin/bash
#SBATCH -J Iso_R1
#SBATCH --nodes 5
#SBATCH --ntasks-per-node=24
#SBATCH --mem=48gb 
#SBATCH --mail-type=all    
#SBATCH -t 03:00:00
#SBATCH -p standard
#SBATCH -o /scratch/dkhaitan/IsolatedS1S2/test_out/matlab.log
ulimit -u 8192
module load matlab
slurm_mdce_start
matlab -nodesktop -r Level1_sbatch >/scratch/dkhaitan/IsolatedS1S2/test_out/log.txt
slurm_mdce_stop
