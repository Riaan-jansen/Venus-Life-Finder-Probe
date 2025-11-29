# Venus-Life-Finder-Probe

2025-11-25

Tools for the analysis of the LEVL probes for the Venus Life Finder mission

## Features

- **Rotation Animation**: Simulate and save animations of rigid body rotations.
- **Multi-View Animation**: Generate animations with multiple perspectives.
- **NDIM root finder for iterative OpenFOAM runs**
- **1D root finder for connect masses rigid body problem**


## Requirements

To run the code, install the required Python packages using the `requirements.txt` file:

```bash
pip install -r requirements.txt
```

### to run the code

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

###

## NDIM root finder
To run:

```bash
python root_finder.py 'OF_data/kateForces.csv'
```

output is to terminal


###

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