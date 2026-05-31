# -*- coding: utf-8 -*-
"""
Created on Wed Apr 10 09:17:41 2024

@author: Kate

riaan changes:

summing contribution of pressure and viscous to forces and moments
"""

import numpy as np
import matplotlib.pyplot as plt
import pandas as pd
import re
import os
import subprocess
import shutil
import sys


# defining functions
# plot forces and moments
def plot_forces_and_moments(file_name):
    case = file_name
    file_path = os.path.join(case, 'postProcessing', 'forces','0','forces.dat')
    
    forceRegex=r"([0-9.Ee\-+]+)\s+\(+([0-9.Ee\-+]+)\s([0-9.Ee\-+]+)\s([0-9.Ee\-+]+)\)\s\(([0-9.Ee\-+]+)\s([0-9.Ee\-+]+)\s([0-9.Ee\-+]+)\)+\s\(+([0-9.Ee\-+]+)\s([0-9.Ee\-+]+)\s([0-9.Ee\-+]+)\)\s\(([0-9.Ee\-+]+)\s([0-9.Ee\-+]+)\s([0-9.Ee\-+]+)\)+"
    t = []
    fpx = []; fpy = []; fpz = []
    fvx = []; fvy = []; fvz = []
    mpx = []; mpy = []; mpz = []
    mvx = []; mvy = []; mvz = []
    pipefile=open(file_path,'r')
    lines = pipefile.readlines()
    for line in lines:
        match=re.search(forceRegex,line)
        if match:
            t.append(float(match.group(1)))
            fpx.append(float(match.group(2)))
            fpy.append(float(match.group(3)))
            fpz.append(float(match.group(4)))
            fvx.append(float(match.group(5)))
            fvy.append(float(match.group(6)))
            fvz.append(float(match.group(7)))
            mpx.append(float(match.group(8)))
            mpy.append(float(match.group(9)))
            mpz.append(float(match.group(10)))
            mvx.append(float(match.group(11)))
            mvy.append(float(match.group(12)))
            mvz.append(float(match.group(13)))

    fig, ax1 = plt.subplots(figsize=(10, 6))

    ax1.set_xlabel('Iteration Number')
    #ax1.set_xlim(300, 2000)
    ax1.set_ylabel('Force (N)')
    #ax1.set_ylim(top = np.max(fpz[400:]))
    ax1.grid()
    ax1.set_title(f'{case} Forces (Pressure)')

    ax1.plot(t, fpx, label='Force x')
    ax1.plot(t, fpy, label='Force y')
    ax1.plot(t, fpz, label='Force z')

    ax2 = ax1.twinx()
    ax2.set_ylabel('Moment (N-m)')
    ax2.plot(t, mpx, linestyle='--', label='Moment x')
    ax2.plot(t, mpy, linestyle='--', label='Moment y')
    ax2.plot(t, mpz, linestyle='--', label='Moment z')
    ax2.set_ylim(top = np.max(mpy))

    fig.tight_layout()
    fig.legend()
    plt.savefig(f'results/{case}_forces_plot.png')


def average_forces(file_name):
    case = file_name
    file_path = os.path.join(case, 'postProcessing', 'forces','0','forces.dat')
    
    forceRegex=r"([0-9.Ee\-+]+)\s+\(+([0-9.Ee\-+]+)\s([0-9.Ee\-+]+)\s([0-9.Ee\-+]+)\)\s\(([0-9.Ee\-+]+)\s([0-9.Ee\-+]+)\s([0-9.Ee\-+]+)\)+\s\(+([0-9.Ee\-+]+)\s([0-9.Ee\-+]+)\s([0-9.Ee\-+]+)\)\s\(([0-9.Ee\-+]+)\s([0-9.Ee\-+]+)\s([0-9.Ee\-+]+)\)+"
    t = []
    fpx = []; fpy = []; fpz = []
    fvx = []; fvy = []; fvz = []
    mpx = []; mpy = []; mpz = []
    mvx = []; mvy = []; mvz = []
    pipefile=open(file_path,'r')
    lines = pipefile.readlines()
    for line in lines:
        match=re.search(forceRegex,line)
        if match:
            t.append(float(match.group(1)))
            fpx.append(float(match.group(2)))
            fpy.append(float(match.group(3)))
            fpz.append(float(match.group(4)))
            fvx.append(float(match.group(5)))
            fvy.append(float(match.group(6)))
            fvz.append(float(match.group(7)))
            mpx.append(float(match.group(8)))
            mpy.append(float(match.group(9)))
            mpz.append(float(match.group(10)))
            mvx.append(float(match.group(11)))
            mvy.append(float(match.group(12)))
            mvz.append(float(match.group(13)))

    if len(fpx) != 2000:
        print(file_name, " is not 2000 iterations long.")
        
    # Assuming fpx, fpy, fpz, fvx, fvy, fvz, mpx, mpy, mpz, mvx, mvy, and mvz are already populated with 2000 float numbers each
    
    # Slice each list to get the last X numbers
    X = 50

    last_500_fpx = fpx[-X:]
    last_500_fpy = fpy[-X:]
    last_500_fpz = fpz[-X:]
    last_500_fvx = fvx[-X:]
    last_500_fvy = fvy[-X:]
    last_500_fvz = fvz[-X:]
    last_500_mpx = mpx[-X:]
    last_500_mpy = mpy[-X:]
    last_500_mpz = mpz[-X:]
    last_500_mvx = mvx[-X:]
    last_500_mvy = mvy[-X:]
    last_500_mvz = mvz[-X:]
    
    # Calculate the average of each list
    average_fpx = sum(last_500_fpx) / len(last_500_fpx)
    average_fpy = sum(last_500_fpy) / len(last_500_fpy)
    average_fpz = sum(last_500_fpz) / len(last_500_fpz)
    average_fvx = sum(last_500_fvx) / len(last_500_fvx)
    average_fvy = sum(last_500_fvy) / len(last_500_fvy)
    average_fvz = sum(last_500_fvz) / len(last_500_fvz)
    average_mpx = sum(last_500_mpx) / len(last_500_mpx)
    average_mpy = sum(last_500_mpy) / len(last_500_mpy)
    average_mpz = sum(last_500_mpz) / len(last_500_mpz)
    average_mvx = sum(last_500_mvx) / len(last_500_mvx)
    average_mvy = sum(last_500_mvy) / len(last_500_mvy)
    average_mvz = sum(last_500_mvz) / len(last_500_mvz)
    
    return average_fpx, average_fpy, average_fpz, average_fvx, average_fvy, average_fvz, average_mpx, average_mpy, average_mpz, average_mvx, average_mvy, average_mvz


def plot_residuals(case):
    # Define the directory path
    directory_path = os.path.join(case, 'logs')
    # Pressure
    file_path = os.path.join(directory_path, 'p_0')

    with open(file_path, 'r') as file:
        data = file.readlines()

    # Process the data
    x_values = []
    p_values = []
    for line in data:
        parts = line.split('\t')
        x_values.append(int(parts[0].strip('s')))
        p_values.append(float(parts[1]))

    # Urel x
    file_path = os.path.join(directory_path, 'Urelx_0')

    with open(file_path, 'r') as file:
        data = file.readlines()

    # Process the data
    x_values = []
    Urelx_values = []
    for line in data:
        parts = line.split('\t')
        x_values.append(int(parts[0].strip('s')))
        Urelx_values.append(float(parts[1]))

    # Urel y
    file_path = os.path.join(directory_path, 'Urely_0')

    with open(file_path, 'r') as file:
        data = file.readlines()

    # Process the data
    x_values = []
    Urely_values = []
    for line in data:
        parts = line.split('\t')
        x_values.append(int(parts[0].strip('s')))
        Urely_values.append(float(parts[1]))

    # Urel z
    file_path = os.path.join(directory_path, 'Urelz_0')

    with open(file_path, 'r') as file:
        data = file.readlines()

    # Process the data
    x_values = []
    Urelz_values = []
    for line in data:
        parts = line.split('\t')
        x_values.append(int(parts[0].strip('s')))
        Urelz_values.append(float(parts[1]))

    # Plot the data
    plt.figure(figsize=(10, 6))
    plt.plot(x_values, p_values, label='p')
    plt.plot(x_values, Urelx_values, label='Urelx')
    plt.plot(x_values, Urely_values, label='Urely')
    plt.plot(x_values, Urelz_values, label='Urelz')

    plt.xlabel('Time (s)')
    plt.ylabel('Value')
    plt.yscale('log')
    plt.title(f'{case} Residuals')
    plt.grid(True)
    plt.legend()

    plt.savefig(f'results/{case}_residuals_plot.png')

# get names of all directories in the folder
files = os.listdir(os.getcwd())
directories = [item for item in files if os.path.isdir(os.path.join(os.getcwd(),item))]

for x in directories:
    if x == 'blank_case':
        directories.remove(x)
    if x == 'meshFolder':
        directories.remove(x)

print("Directories in the current directory (excluding 'blank_case'):", directories)

os.mkdir('results')
 
pitch_angles = []
cone_angles = []
omegas = []

fx = []; fy = []; fz = []
mx = []; my = []; mz = []

for name in directories:
    try:
        variables = name.split('_')
        cone_angles.append(variables[0])
        omegas.append(variables[1])
        pitch_angles.append(variables[2])

        # make plots of forces, moments, residuals and save under case name
        plot_residuals(name)
        plot_forces_and_moments(name)
    
        average_fpx, average_fpy, average_fpz, average_fvx, average_fvy, average_fvz, average_mpx, average_mpy, average_mpz, average_mvx, average_mvy, average_mvz = average_forces(name)

        fx.append(np.array(average_fpx) + np.array(average_fvx))
        fy.append(np.array(average_fpy) + np.array(average_fvy))
        fz.append(np.array(average_fpz) + np.array(average_fvz))

        mx.append(np.array(average_mpx) + np.array(average_mvx))
        my.append(np.array(average_mpy) + np.array(average_mvy))
        mz.append(np.array(average_mpz) + np.array(average_mvz))

    except Exception as e:
        print(f'Error: {e}\n(Directory name: "{name}" does not fit pattern)')
        directories.remove(name)

# get vz
parent_folder = os.path.basename(os.getcwd())

parent_name = parent_folder.split('_')

vz_flt = float(parent_name[-1])

vz = np.ones(len(directories))*vz_flt

print("len dir", len(directories))
print(vz[0], vz[-1])
print(len(mx))

# Create DataFrame
# comment on case so header doesnt get read in as numeric
data = {
    '# Case': directories,
    'Cone Angle': cone_angles,
    'Omega': omegas,
    'Pitch Angle': pitch_angles,
    'vz': vz,
    'Average Fx': fx,
    'Average Fy': fy,
    'Average Fz': fz,
    #'Average Fvx': fvx,
    #'Average Fvy': fvy,
    #'Average Fvz': fvz,
    'Average Mx': mx,
    'Average My': my,
    'Average Mz': mz
    #'Average Mvx': mvx,
    #'Average Mvy': mvy,
    #'Average Mvz': mvz
}

df = pd.DataFrame(data)

# Save DataFrame to CSV file in results directory
df.to_csv('results/forces_results.csv', index=False)
