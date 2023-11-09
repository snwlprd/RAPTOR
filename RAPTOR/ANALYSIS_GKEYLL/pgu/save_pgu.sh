#!/bin/bash
#SBATCH --job-name=energies
#SBATCH --ntasks=16
#SBATCH --time=01-00:00:00
#SBATCH --mem-per-cpu=3000M
#SBATCH --output=/nfs/scratch/edyveaja/GYKL_RUNS/sims/MR_2048/analysis/pgu.log
#SBATCH --partition=gpu

source activate Mass_Ratios
python3 pi_d_pTh_Eth_pgu.py 



