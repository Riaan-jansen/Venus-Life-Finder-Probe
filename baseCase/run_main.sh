#!/bin/bash

#SBATCH -p stellarprop             # Specifies that the job will run on the default queue nodes.
#SBATCH --job-name=main          # A name for the job to be used when Job Monitoring
#SBATCH -D .                            # set wd
#SBATCH --time=168:00:00                    # maximum run time for a job hrs:min:sec
#SBATCH --nodes=1                           # Number of full nodes to use
#SBATCH --ntasks-per-node=1                 # Run a single task on each node
#SBATCH --ntasks=1                          # Run a single task (increase this value for parallelisation across CPUs)
#SBATCH -A Research_Project-129871
#SBATCH --output=test-single-core-%j.out    # Standard output and error log
#SBATCH --mail-type=END                     # Mail events (NONE, BEGIN, END, FAIL, ALL)
#SBATCH --mail-user=rj429@exeter.ac.uk      # E-Mail address of the user that needs to be notified.

# Load Modules
module load GCC/12.3.0
module load OpenMPI/4.1.5-GCC-12.3.0

# for postprocessing
module use /lustre/shared/isca_compute/modules/all
module load matplotlib/3.7.0-gfbf-2022b

# Load OpenFOAM Environment
# We source the bashrc file AFTER we load the modules
# source /lustre/home/rj429/OpenFOAM/OpenFOAM-10/etc/bashrc

# run
python main.py
