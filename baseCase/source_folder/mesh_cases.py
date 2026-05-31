# -*- coding: utf-8 -*-
"""
Created on Thu Mar 21 12:37:30 2024

@author: Kate

riaan changes:

now the script is executed as part of the main.py script, takes parameters from
entries in a "parameters.txt" file which are passed as arguments when calling
this script.

to run individually (but you probably dont want to do that):

 - python make_and_run_cases.py "{cone}" "{omegas}" pitch

to run everything (recommended):

 - python main.py

"""

import os
import subprocess
import shutil
import sys
import json


# parameters as read from parameters.txt
cone_angles = sys.argv[1]
omegas = sys.argv[2]  		# Measured in RPM
pitch_angle = sys.argv[3]
vz = sys.argv[4]                # in m/s

# convert format
cone_angles = json.loads(cone_angles)
omegas = json.loads(omegas)
pitch_angle = json.loads(pitch_angle)
vz = json.loads(vz)

cone_angles = [int(x) for x in cone_angles]
omegas = [int(x) for x in omegas]
pitch_angle = [int(x) for x in pitch_angle]
vz = float(vz)

# variable pass check
print(f"""############# MESHING FOR #############
 pitch = {pitch_angle}\n cone = {cone_angles}\n
######################################\n
""")

# specify the folder to store meshes in (AND must contain original_blade.stl)
blank_folder = 'blank_case' # "copy_case" # 2024-05-04

# Checks to see that this case exists in the directory that the python script is running from
if os.path.isdir(blank_folder) is False:
    print("Error: Path ", blank_folder, " does not exist. Exiting script.")
    sys.exit(0)
else:
    pass

print("Pitch angle = ", pitch_angle)

# loop through cone angle and omegas
for cone in cone_angles: # 2024-05-04
    for pitch in pitch_angle: # 2024-05-04
        print("Meshing for case: Cone angle = ", cone, " pitch angle = ", pitch)

        src = blank_folder
        dst = f"geom_{cone}_{pitch}"

        # check that copied case does NOT already exist (to ensure a simulation isn't overwritten)
        if os.path.isdir(dst):
            print("Error: Path ", dst, " already exists. Remove the directory or rename it to run this script.")
            continue # skips to next angle if this angle exists
        # copy the original case
        else:

############################ CASE FOLDER MAKING ####################################

            try:
                shutil.copytree(src, dst)
                print(f"Directory copied successfully from '{src}' to '{dst}'.")
            except Exception as e:
                print(f"Error in copying source folder: {e}")
                sys.exit(0)
                _ = input('foo') # continue 2024-05-04


############################ initialConditions FILE MAKING ####################################
                # need to make here first so decomposePar doesnt have a fit

                print("Making initialConditions file")
                initConds_contents = fr"""
/*--------------------------------*- C++ -*----------------------------------*\
=========                 |
\\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
 \\    /   O peration     | Website:  https://openfoam.org
  \\  /    A nd           | Version:  10
   \\/     M anipulation  |
\*---------------------------------------------------------------------------*/

flowVelocity         (0 0 {vz});
pressure             0;
turbulentKE          0.24;
turbulentOmega       1.78;
magU                    {vz};
// ************************************************************************* //
"""

                init_path = os.path.join(dst,"0/include","initialConditions")
                with open(init_path, 'w') as file:
                    file.write(initConds_contents)
                os.chmod(init_path,0o644)

                if os.path.exists(init_path) is False:
                    print("Error: Path ", init_path, " does not exist. Exiting script.")
                    sys.exit(0)
                else: pass


############################ ALLRUN FILE MAKING ####################################

            print("Making Allrun file")
            allrun_file_contents = r"""#!/bin/sh
set -e
cd ${0%/*} || exit 1    # Run from this directory

# Source tutorial run functions
# . $WM_PROJECT_DIR/bin/tools/RunFunctions

# apparently necessary for mpi on ISCA
export OMPI_MCA_btl='^uct,ofi'
export OMPI_MCA_pml='ucx'
export OMPI_MCA_mtl='^ofi'
"""

            allrun_file_contents +="""# rotating the probe
surfaceTransformPoints -rotate-x {1} geometry/Simplified_Probe_smooth.stl constant/triSurface/blade_1.stl
surfaceTransformPoints -rotate-y -{0} constant/triSurface/blade_1.stl constant/triSurface/blade.stl

surfaceFeatureExtract
blockMesh

snappyHexMesh -overwrite
checkMesh -allGeometry

#------------------------------------------------------------------------------""".format(cone, pitch)
            allrun_path = os.path.join(dst,"Allrun")
            with open(allrun_path, 'w') as file:
                file.write(allrun_file_contents)

            if os.path.exists(allrun_path) is False:
                print("Error: Path ", allrun_path, " does not exist. Exiting script.")
                sys.exit(0)
            else:
                pass

            # modify the rights of the file, to give the user run rights
            os.chmod(allrun_path,0o755)

            print("successfully made allrun file")

            # Run the allrun bash script
            Allrun = dst+f"/Allrun"# > {dst}/log_allrun &"
            subprocess.call(Allrun, shell=True)

print("Completed meshing - moving meshes to folder")

os.system("mkdir meshFolder && mv geom_* meshFolder/.")

