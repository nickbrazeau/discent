#!/bin/bash
#SBATCH --job-name=polysimIBD      # Job name
#SBATCH --mail-type=END,FAIL         # Mail events (NONE, BEGIN, END, FAIL, ALL)
#SBATCH --mail-user=nbrazeau@med.unc.edu
#SBATCH --ntasks=24                  # Number of MPI tasks (i.e. processes)
#SBATCH --mem-per-cpu=600mb          # Memory (i.e. RAM) per processor
#SBATCH --time=36:00:00              # Wall time limit (days-hrs:min:sec)
#SBATCH --output=polysim_%j.log     # Path to the standard output and error files

R CMD BATCH run.R

# Removing .RData is recommended.
rm -f .RData
