#!/bin/bash
#SBATCH --job-name=polysimIBD      # Job name
#SBATCH --mail-type=END,FAIL         # Mail events (NONE, BEGIN, END, FAIL, ALL)
#SBATCH --mail-user=nbrazeau@med.unc.edu
#SBATCH --nodes=1                    # Run all processes on a single node
#SBATCH --ntasks=1                  # Number of MPI tasks (i.e. processes)
#SBATCH --cpus-per-task=24
#SBATCH --mem=64G          # Memory (i.e. RAM) per processor
#SBATCH --time=36:00:00              # Wall time limit (days-hrs:min:sec)
#SBATCH --output=polysim_%j.log     # Path to the standard output and error files

R CMD BATCH _future_polySimIBD.R

