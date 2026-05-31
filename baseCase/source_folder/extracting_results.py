# -*- coding: utf-8 -*-
"""
Created on Wed Apr 10 09:17:41 2024

@author: Kate
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
def read_lines(case_name, file):
    """input: case directory name and file (force.dat or moment.dat)
    returns: t, fpx, fpy, fpz, fvx, fvy, fvz"""
    case = case_name
    # takes you to the forces given called from the case directory ie the file
    # of interest (sim results) is in the same directory you are calling from
    file_path = os.path.join(case, 'postProcessing', 'forces', '0', file)  # 10-10-25 forces.dat -> force.dat

    # forceRegex=r"([0-9.Ee\-+]+)\s+\(+([0-9.Ee\-+]+)\s([0-9.Ee\-+]+)\s([0-9.Ee\-+]+)\)\s\(([0-9.Ee\-+]+)\s([0-9.Ee\-+]+)\s([0-9.Ee\-+]+)\)+\s\(+([0-9.Ee\-+]+)\s([0-9.Ee\-+]+)\s([0-9.Ee\-+]+)\)\s\(([0-9.Ee\-+]+)\s([0-9.Ee\-+]+)\s([0-9.Ee\-+]+)\)+"
    # forceRegex=r"^\d+\s+(?:[-+]?\d+\.\d+e[-+]?\d+\s){9}$"
    t = []
    fpx = []; fpy = []; fpz = []
    fvx = []; fvy = []; fvz = []

    pipefile=open(file_path,'r')
    lines = pipefile.readlines()
    # 10-10-25 replaced regex searching which wasnt working with split lines.
    # relies on beginning of data being line 4! not ideal
    

    matches = 0
    for line in lines:

        #match=re.search(forceRegex,line)
        matches = matches + 1
        line = line.split()
        if line[0] != '#':
            t.append(float(line[0]))
            fpx.append(float(line[1]))
            fpy.append(float(line[2]))
            fpz.append(float(line[3]))
            # fvx.append(float(line[4]))
            # fvy.append(float(line[5]))
            # fvz.append(float(line[6]))
            fvx.append(float(line[7]))
            fvy.append(float(line[8]))
            fvz.append(float(line[9]))
            # mvx.append(float(line[10]))
            # mvy.append(float(line[11]))
            # mvz.append(float(line[12]))

    return t, fpx, fpy, fpz, fvx, fvy, fvz


# plot forces and moments
def plot_forces_and_moments(case_name):

    t, fpx, fpy, fpz, fvx, fvy, fvz = read_lines(case_name, 'force.dat')
    t, mpx, mpy, mpz, mvx, mvy, mvz = read_lines(case_name, 'moment.dat')
    
    fig, ax1 = plt.subplots(figsize=(10, 6))

    ax1.set_xlabel('Iteration Number')
    # ax1.set_xlim(0, 800)
    ax1.set_ylabel('Force (N)')
    #ax1.set_ylim(top = np.max(fpz[400:]))
    ax1.grid()
    ax1.set_title(f'{case_name} Forces (Pressure)')

    x1 = 2000
    t = t[:x1]
    fpx = fpx[:x1]
    fpy = fpy[:x1]
    fpz = fpz[:x1]
    mpx = mpx[:x1]
    mpy = mpy[:x1]
    mpz = mpz[:x1]

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
    plt.savefig(f'results/{case_name}_forces_plot.png')


def average_forces(case_name):

    t, fpx, fpy, fpz, fvx, fvy, fvz = read_lines(case_name, 'force.dat')
    t, mpx, mpy, mpz, mvx, mvy, mvz = read_lines(case_name, 'moment.dat')
    

    if len(fpx) != 2000:
        print(case_name, " is not 2000 iterations long.")
        
    # Assuming fpx, fpy, fpz, fvx, fvy, fvz, mpx, mpy, mpz, mvx, mvy, and mvz are already populated with 2000 float numbers each
    
    # Slice each list to get the last 500 numbers
    # 10-10-25  This only works if results reached convergence in the last 500 iterations
    last_500_fpx = fpx[-500:]
    last_500_fpy = fpy[-500:]
    last_500_fpz = fpz[-500:]
    last_500_fvx = fvx[-500:]
    last_500_fvy = fvy[-500:]
    last_500_fvz = fvz[-500:]
    last_500_mpx = mpx[-500:]
    last_500_mpy = mpy[-500:]
    last_500_mpz = mpz[-500:]
    last_500_mvx = mvx[-500:]
    last_500_mvy = mvy[-500:]
    last_500_mvz = mvz[-500:]
    
    # Calculate the average of each list
    average_fpx = sum(last_500_fpx) / len(last_500_fpx)  # is this not just 500?
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
    plt.close()

# get names of all directories in a file:

files = os.listdir(os.getcwd())
directories = [item for item in files if os.path.isdir(os.path.join(os.getcwd(),item)) and item != 'blank_case' and item != 'meshFolder']
# directories = ['10_240_-5']
# added "item != blank_case" to not loop twice also avoids string literal error 3/4/25
# for 'blank_case' in directories:
#     directories.remove('blank_case')
# also changed copy_case to blank_case
    
print("Directories in the current directory (excluding 'blank_case'):", directories)

os.mkdir('results')

#df_results = pd.DataFrame(directories, columns=['Names'])

#for index, row in df_results.iterrows():
    
    # separate name into cone angle, pitch angle, rotation rate
 
pitch_angles = []
cone_angles = []
omegas = []
fpx = []; fpy = []; fpz = []
fvx = []; fvy = []; fvz = []
mpx = []; mpy = []; mpz = []
mvx = []; mvy = []; mvz = []


for name in directories:
    
    variables = name.split('_')
    cone_angles.append(variables[0])
    omegas.append(variables[1])
    pitch_angles.append(variables[2])
    
    #plot_residuals(name)
    #plot_forces_and_moments(name)
    
    average_fpx, average_fpy, average_fpz, average_fvx, average_fvy, average_fvz, average_mpx, average_mpy, average_mpz, average_mvx, average_mvy, average_mvz = average_forces(name)
    
    fpx.append(average_fpx)
    fpy.append(average_fpy)
    fpz.append(average_fpz)
    mpx.append(average_mpx)
    mpy.append(average_mpy)
    mpz.append(average_mpz)


# Create DataFrame
data = {
    'Case': directories,
    'Cone Angle': cone_angles,
    'Omega': omegas,
    'Pitch Angle': pitch_angles,
    'Average Fx': fpx,
    'Average Fy': fpy,
    'Average Fz': fpz,
    'Average Mx': mpx,
    'Average My': mpy,
    'Average Mz': mpz,
}

df = pd.DataFrame(data)

print("df to csv okay")
# Save DataFrame to CSV file in results directory
df.to_csv('results/forces_results.csv', index=False)
    
    
for name in directories:
    plot_residuals(name)
    plot_forces_and_moments(name)
