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

to run everything:

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
pitch_angle = int(pitch_angle)
vz = float(vz)

# variable pass check
print(f"""############# PARAMETER CHECK #############
 pitch = {pitch_angle}\n cone = {cone_angles}\n omegas = {omegas}.
######################################\n
""")

# specify the 'base case' you wish to use
original_case = 'blank_case' # "copy_case" # 2024-05-04

# Checks to see that this case exists in the directory that the python script is running from
if os.path.isdir(original_case) is False:
    print("Error: Path ", original_case, " does not exist. Exiting script.")
    sys.exit(0)
else: pass

print("Pitch angle = ", pitch_angle)

# loop through cone angle and omegas
for angle in cone_angles: # 2024-05-04
    for omega in omegas: # 2024-05-04
        print("Cone angle = ", angle)

        # copy original_case as cone_angle_30 e.g.
        # from https://pynative.com/python-copy-files-and-directories/#h-copy-entire-directory :
        # shutil.copytree(src, dst, symlinks=False, ignore=None, copy_function=copy2, ignore_dangling_symlinks=False, dirs_exist_ok=False
        # define source src and destination dst
        # The dirs_exist_ok dictates whether to raise an exception in case dst or any missing parent directory already exists.
        src = original_case
        dst = "{}_{}_{}".format(angle,omega,pitch_angle)
        
        # check that copied case does NOT already exist (to ensure a simulation isn't overwritten)
        if os.path.isdir(dst):
            print("Error: Path ", dst, " already exists. Remove the directory or rename it to run this script.")
            #continue # skips to next angle if this angle exists
        # copy the original case
        else:

############################ CASE FOLDER MAKING ####################################

            try:
                shutil.copytree(src, dst)
                print(f"Directory copied successfully from '{src}' to '{dst}'.")
            except Exception as e:
                print(f"An error occurred: {e}")
                sys.exit(0)
                _ = input('foo') # continue 2024-05-04


############################ SRFProperties file MAKING ####################################

            print("Making SRFProperties File")
            srf_file_contents = r"""
            /*--------------------------------*- C++ -*----------------------------------*\
              =========                 |
              \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
               \\    /   O peration     | Website:  https://openfoam.org
                \\  /    A nd           | Version:  10
                 \\/     M anipulation  |
            \*---------------------------------------------------------------------------*/
            FoamFile
            {
                format      ascii;
                class       dictionary;
                location    "constant";
                object      SRFProperties;
            }
            // * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //
            
            SRFModel        rpm;
            
            origin          (0 0 -5);
            axis            (0 0 -1);
            
            rpmCoeffs
            {
                rpm """            
            srf_file_contents +="{}".format(omega)
            srf_file_contents +=""";
            }
            
            
            // ************************************************************************* //"""

            srf_path = os.path.join(dst,"constant","SRFProperties")
            with open(srf_path, 'w') as file:
                file.write(srf_file_contents)
            # copied from the motorbike tutorial but edited
            #file.close()
            os.chmod(srf_path,0o644)
            
            if os.path.exists(srf_path) is False:
                print("Error: Path ", srf_path, " does not exist. Exiting script.")
                sys.exit(0)
            else: pass


############################ initialConditions FILE MAKING ####################################

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
cd ${0%/*} || exit 1    # Run from this directory

# Source tutorial run functions

. $WM_PROJECT_DIR/bin/tools/RunFunctions

# apparently necessary for mpi on ISCA
export OMPI_MCA_btl='^uct,ofi'
export OMPI_MCA_pml='ucx'
export OMPI_MCA_mtl='^ofi'
"""

            allrun_file_contents +="""# rotating the probe
runApplication surfaceTransformPoints "Rx={1}" constant/triSurface/original_blade.stl constant/triSurface/blade_step1.stl
mv log.surfaceTransformPoints log.surfaceTransformPoints1
runApplication surfaceTransformPoints "Ry={0}" constant/triSurface/blade_step1.stl constant/triSurface/blade.stl

# extracting the surface features *.eMesh + create a boundary
runApplication surfaceFeatures
runApplication blockMesh

# parallel meshing
runApplication decomposePar
runParallel snappyHexMesh -overwrite
runApplication reconstructParMesh -constant

# decomposePar again for the fields
rm -rf processor*
runApplication checkMesh
mv log.decomposePar log.decomposePar1
runApplication decomposePar

# run simplefoam
runParallel SRFSimpleFoam
### runParallel SRFSimpleFoam -postProcess -func yPlus -latestTime

# reconstruct
runApplication reconstructParMesh -constant
runApplication reconstructPar -latestTime

runApplication foamLog log.SRFSimpleFoam
rm -r processor*
#------------------------------------------------------------------------------""".format(angle, pitch_angle)
            allrun_path = os.path.join(dst,"Allrun")
            with open(allrun_path, 'w') as file:
                file.write(allrun_file_contents)

            if os.path.exists(allrun_path) is False:
                print("Error: Path ", allrun_path, " does not exist. Exiting script.")
                sys.exit(0)
            else: pass

		# modify the rights of the file, to give the user run rights
            os.chmod(allrun_path,0o755)

            print("successfully made allrun file")

		# Run the allrun bash script
            Allrun = dst+f"/Allrun"# > {dst}/log_allrun &"
            subprocess.call(Allrun, shell=True)
            print(f"Ran Allrun for a cone angle ({angle}) degrees, omega ({omega}), pitch angle ({pitch_angle}). Check log for more details")

print("Completed script for pitch angle = ", pitch_angle)
