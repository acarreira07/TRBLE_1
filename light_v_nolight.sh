#!/bin/bash
#SBATCH --job-name=light_v_nolight.
#SBATCH --time=4-00:00:00
#SBATCH --ntasks=1
#SBATCH --partition=trc
#SBATCH --output=./logs/mainlog.out
#SBATCH --open-mode=append
#SBATCH --cpus-per-task=4
#SBATCH --mail-type=ALL

ml R
R -u /home/users/arnaldo/TRBLE_1/treble_larvae_behavior_walkthrough_mz1407.R 
