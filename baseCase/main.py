# -*- coding: utf-8 -*-
"""
Created on Sat Dec 20 19:37:30 2025

@author: riaan

To run:

 - python main.py

What it does:

 - runs openfoam simulations for all of the variables in "parameters.txt"
 - runs for every pitch angle and (will) process all of the results
    - for now results can be processed from the base directory with cat_results.sh
"""

import sys
import shutil
import os
import subprocess

# number of cores to use, NB: also needs to be changed in source_folder/blank_case/system/decomposeParDict
nCores = 16

# reading the parameters file
with open("parameters.txt", "r") as f:
    parameters = f.readlines()

count = 0
for line in parameters:
    # comment line - leave no non-data line uncommented
    if line[0] == '#':
        pass

    else:
        count = count + 1

        # means that data has to be in a specific order
        if count == 1:
            pitch = line.split(',')
            pitch = [float(x) for x in pitch]

        elif count == 2:
            omega = line.split(',')
            omega = [float(x) for x in omega]

        elif count == 3:
            cone = line.split(',')
            cone = [float(x) for x in cone]

        elif count == 4:
            vz = line.split(',')
            vz = [float(x) for x in vz]


# base directory main is called from and source_folder is copied from
base_dir = os.path.dirname(os.path.abspath(__file__))

# makes a directory and executes a job for every pitch angle
for v in vz:

    # change to base directory
    os.chdir(base_dir)

    src = "source_folder"
    dst = f"vz_{int(v*100)}"

    # check that copied case does NOT already exist (to ensure a simulation isn't overwritten)
    if os.path.isdir(dst):
        print("Error: Path ", dst, " already exists. Remove the directory or rename it to run this script.")

    # copy the original case
    else:
        try:
            shutil.copytree(src, dst)
            print(f"Directory copied successfully from '{src}' to '{dst}'.")
        except Exception as e:
            print(f"An error occurred: {e}")
            sys.exit(0)
            _ = input('foo') # continue 2024-05-04

    # change to select pitch directory
    os.chdir(dst)

    # writing the slurm job
    job_text=f"""#!/bin/bash

#SBATCH -p stellarprop             # Specifies that the job will run on the default queue nodes.
#SBATCH --job-name=vz{int(v*100)}          # A name for the job to be used when Job Monitoring
#SBATCH -D .                            # set wd
#SBATCH --time=168:00:00                    # maximum run time for a job hrs:min:sec
#SBATCH --nodes={int(nCores/16)}                           # Number of full nodes to use
#SBATCH --ntasks-per-node=16                 # Run a single task on each node
#SBATCH --ntasks={nCores}                          # Run a single task (increase this value for parallelisation across CPUs)
#SBATCH -A Research_Project-129871
#SBATCH --output=test-single-core-%j.out    # Standard output and error log
#SBATCH --mail-type=END                     # Mail events (NONE, BEGIN, END, FAIL, ALL)
#SBATCH --mail-user=rj429@exeter.ac.uk      # E-Mail address of the user that needs to be notified.

# Load Modules
module load GCC/12.3.0
module load OpenMPI/4.1.5-GCC-12.3.0
module load matplotlib/3.7.0-gfbf-2022b

# for postprocessing
module use /lustre/shared/isca_compute/modules/all

# Load OpenFOAM Environment
# We source the bashrc file AFTER we load the modules
# source /lustre/home/rj429/OpenFOAM/OpenFOAM-10/etc/bashrc
module load OpenFOAM/v2406-foss-2023a

source $FOAM_BASH

# mesh
python mesh_cases.py "{cone}" "{omega}" "{pitch}" "{v}"

# run
python run_cases.py "{cone}" "{omega}" "{pitch}" "{v}"

# results
python extracting_results.py
"""

    # writing the job script
    # job_path = os.path.join(dst,"job.sh")
    job_path = 'job.sh'
    with open(job_path, 'w') as file:
        file.write(job_text)

        if os.path.exists(job_path) is False:
            print("Error: Path ", allrun_path, " does not exist. Exiting script.")
            sys.exit(0)
        else: pass

        # modify the rights of the file, to give the user run rights
        os.chmod(job_path,0o755)

    # Run the allrun bash script
    Allrun = "sbatch job.sh"
    subprocess.call(Allrun, shell=True)

print("All pitch angles began running")
