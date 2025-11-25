# -*- coding: utf-8 -*-
"""
Created on Thu Mar 21 12:37:30 2024

@author: Kate

NJ comments
    2024-05-04
        KS using openfoam10
        Now installed on Acer01
        On booting Ubuntu, type 'openfoam10' to "source the environment"
            ie set environment variables. Type 'env' to see the result
        
        original_blade does not exist - copied from blade.stl (which appears to be at a jaunty angle)
        
        Both surfaceTransformPoints calls rotate about Ry - is this intended?
        
        snappyHexMesh crashes: decomposePar moved after checkMesh
        
        blockMesh size reduced by 4x in linear directions, # cells reduced by 4**3 = 64
            If this doesn't process...
        TODO: put wing in more sensible location in volume - one third or one quarter "in"
        
        TODO: Have had to "tone down" mesh considerably - check KS & my results similar
        TODO: If changing maxLocalCells & ditto Global lets sHM run, try changing other mods back
"""


import os
import subprocess
import shutil
import sys

# specify what parameters you want investigated. Only 1 pitch angle per script can be ran, but as many cone angles and omegas as you wish

pitch_angle = -5 
cone_angles = [0, 10, 5]
omegas = [480, 240, 360] # Measured in RPM
vz = 2.5 # purely for printing

# specify the 'base case' you wish to use
original_case = 'blank_case' # "copy_case" # 2024-05-04


# Checks to see that this case exists in the directory that the python script is running from
if os.path.isdir(original_case) is False:
    print("Error: Path ", original_case, " does not exist. Exiting script.")
    sys.exit(0)
else: pass

print("Pitch angle = ",pitch_angle)

# loop through cone angle and omegas
# for angle in cone_angles: # 2024-05-04
angle = cone_angles[1]
if 1:
    # for omega in omegas: # 2024-05-04
    omega = omegas[1]
    if 1:
        print("Cone angle = ", angle)
        ratio = vz/omega
        print("vz/omega ratio =", vz/omega)
        # copy original_case as cone_angle_30 e.g.
        # from https://pynative.com/python-copy-files-and-directories/#h-copy-entire-directory :
        # shutil.copytree(src, dst, symlinks=False, ignore=None, copy_function=copy2, ignore_dangling_symlinks=False, dirs_exist_ok=False
        # define source src and destination dst
        # The dirs_exist_ok dictates whether to raise an exception in case dst or any missing parent directory already exists.
        src = original_case
        dst = f"{angle}_{omega}_{pitch_angle}"

        # check that copied case does NOT already exist (to ensure a simulation isn't overwritten)
        if os.path.isdir(dst):
            print("Error: Path ", dst, " already exists. Remove the directory or rename it to run this script.")
            #continue # skips to next angle if this angle exists
        else:
            # copy the original case
            try:
                shutil.copytree(src, dst)
                print(f"Directory copied successfully from '{src}' to '{dst}'.")
            except Exception as e:
                print(f"An error occurred: {e}")
                sys.exit(0)

            print("Making SRFProperties File")
            srf_file_contents = fr"""
            /*--------------------------------*- C++ -*----------------------------------*\
              =========                 |
              \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
               \\    /   O peration     | Website:  https://openfoam.org
                \\  /    A nd           | Version:  10
                 \\/     M anipulation  |
            \*---------------------------------------------------------------------------*/
            FoamFile
            {{
                format      ascii;
                class       dictionary;
                location    "constant";
                object      SRFProperties;
            }}
            // * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //
            
            SRFModel        rpm;
            
            origin          (0 0 -5);
            axis            (0 0 -1);
            
            rpmCoeffs
            {{
                rpm {omega};
            }}
            
            
            // ************************************************************************* //"""

            # writing the SRFProperties file
            srf_path = os.path.join(dst,"constant","SRFProperties")
            with open(srf_path, 'w') as file:
                file.write(srf_file_contents)
            # copied from the motorbike tutorial but edited
            os.chmod(srf_path,0o644)
            
            if os.path.exists(srf_path) is False:
                print("Error: Path ", srf_path, " does not exist. Exiting script.")
                sys.exit(0)
            else: pass
            
            # making the Allrun file
            print("Making Allrun file")
            allrun_file_contents = fr"""#!/bin/sh
            cd ${{0%/*}} || exit 1    # Run from this directory
            
            # Source tutorial run functions
            . $WM_PROJECT_DIR/bin/tools/RunFunctions
            
	    surfaceTransformPoints -rotate-y {pitch_angle} constant/triSurface/original_blade.stl constant/triSurface/blade_1.stl
	    mv log.surfaceTransformPoints log.surfaceTransformPoints1
	    surfaceTransformPoints -rotate-x {angle} constant/triSurface/blade_1.stl constant/triSurface/blade.stl
            surfaceFeatureExtract
	    blockMesh

	    snappyHexMesh -overwrite
	    checkMesh -allGeometry
	    decomposePar
	    mpirun -np 16 SRFSimpleFoam -parallel -postProcess -func yPlus > log.SRFSimpleFoam
	    reconstructParMesh -constant
	    reconstructPar -latestTime

	    foamLog log.SRFSimpleFoam
	    rm -r processor*
            #------------------------------------------------------------------------------"""

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
            print(f"Ran Allrun for a cone angle ({angle}) degrees, omega ({omega}), pitch angle ({pitch_angle}), velocity/omega ratio ({ratio}). Check log for more details")

print("Completed script for pitch angle = ", pitch_angle)
