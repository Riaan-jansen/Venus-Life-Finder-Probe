# Venus-Life-Finder-Probe

2025-11-25

Tools for the analysis of the LEVL probes for the Venus Life Finder mission

## Folders

- **MINIMISER**: contains minimisation script and plotting tools
- **cfd_data**: data collated in a format that works with the analysis tools
- **baseCase**: final working openfoam simulation settings and running script
- **other**: anything else is not important

## Requirements

required to run minimiser.py:

- import numpy as np
- from scipy.interpolate import LinearNDInterpolator
- from scipy.optimize import least_squares
- from scipy.spatial.transform import Rotation as R
- import matplotlib.pyplot as plt
- import argparse


### To run the minimal solution finder

```bash
python minimiser.py <filepath> <scaling> <window>
```

where filepath is the path to a valid .csv file containing openfoam results
and scaling is the weighting applied to "force" over moments (this is to 
compare similar magnitudes).
Window argument tells the program what simulations to take into account
(ignoring those with insufficient lift).

# output

output is to terminal (can be piped into file with $ ... > log )

firstly reads the minimal combination from the table of ran simulations

then reads the interpolated solution and associated residual forces and moments
from net zero (how far away from steady state).

### miscellaneous

## rotating probe animation

# first run the rigid body solver

```bash
gfortran -o <exe_name> omega_solver.f90
./<exe_name> > <data_file>
```

# run python plotting
python ani_omega.py <data_file> [plot_option]

plot options are 0: axes attached single view; 1: multi-view
saves in /results/.


## 1D root finder to two masses on a rod problem

```bash
gfortran -o trig trig_solver.f90
./trig
```

automatically writes to 'datadump.txt'

to plot:

```bash
python plot_probe.py datadump.txt
```