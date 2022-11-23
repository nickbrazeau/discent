#!/bin/bash
#SBATCH --job-name=polysimIBD      # Job name
#SBATCH --mail-type=END,FAIL         # Mail events (NONE, BEGIN, END, FAIL, ALL)
#SBATCH --mail-user=nbrazeau@med.unc.edu
#SBATCH --ntasks=24                  # Number of MPI tasks (i.e. processes)
#SBATCH --cpus-per-task=1            # Number of cores per MPI task
#SBATCH --nodes=1028                    # Maximum number of nodes to be allocated
#SBATCH --ntasks-per-node=12         # Maximum number of tasks on each node
#SBATCH --ntasks-per-socket=6        # Maximum number of tasks on each socket
#SBATCH --distribution=cyclic:cyclic # Distribute tasks cyclically first among nodes and then among sockets within a node
#SBATCH --mem-per-cpu=600mb          # Memory (i.e. RAM) per processor
#SBATCH --time=36:00:00              # Wall time limit (days-hrs:min:sec)
#SBATCH --output=polysim_%j.log     # Path to the standard output and error files

R CMD BATCH run.R

# Removing .RData is recommended.
rm -f .RData
