#!/bin/bash
# Array job: one task = one simulated table fitted in all six parameterisations (see R/sim_tables.R).
# 3 prior widths x 3 penalties x 30 tables = 270 tasks. Edit the lines marked CHECK for Hamilton.
#SBATCH --job-name=sim_tables
#SBATCH --array=1-270
#SBATCH --cpus-per-task=2
#SBATCH --mem=4G
#SBATCH --time=02:00:00
#SBATCH --output=logs/sim_%A_%a.out
##SBATCH --partition=...        # CHECK: partition / account lines for Hamilton, if needed

# CHECK: however you load R on Hamilton, e.g.  module load R   (use the same R and cmdstanr as in RStudio)
mkdir -p logs
cd "$SLURM_SUBMIT_DIR"           # submit from the repository root
Rscript scripts/run_sim_task.R
