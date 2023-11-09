#!/bin/bash
#SBATCH --job-name=energies
#SBATCH --ntasks=16
#SBATCH --time=05:00:00
#SBATCH --mem-per-cpu=3000M
#SBATCH --output=/nfs/scratch/edyveaja/GYKL_RUNS/sims/MR_2048/analysis/en.log


source activate Mass_Ratios
python3 lpg_save_energies.py



