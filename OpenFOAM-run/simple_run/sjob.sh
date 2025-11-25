#!/bin/bash

#SBATCH -p stellarprop             # Specifies that the job will run on the default queue nodes.
#SBATCH --job-name=paraprobe          # A name for the job to be used when Job Monitoring 
#SBATCH -D . 				# set wd
#SBATCH --time=168:00:00                    # maximum run time for a job hrs:min:sec
#SBATCH --nodes=1                           # Number of full nodes to use
#SBATCH --ntasks-per-node=16                 # Run a single task on each node
#SBATCH --ntasks=16                          # Run a single task (increase this value for parallelisation across CPUs)
#SBATCH -A Research_Project-129871
#SBATCH --output=test-single-core-%j.out    # Standard output and error log
#SBATCH --mail-type=END                     # Mail events (NONE, BEGIN, END, FAIL, ALL)
#SBATCH --mail-user=rj429@exeter.ac.uk      # E-Mail address of the user that needs to be notified.


# create and set appropriate temp dir variable for your software
# you will need to uncomment the appropriate line, adjust the actual directory name
# and read your software manual to figure out which variable is appropriate.

#load modules
module use /lustre/shared/isca_compute/modules/all
module load OpenFOAM/v2406-foss-2023a

source $FOAM_BASH


#run
#./allRun
python make_and_run_cases.py
