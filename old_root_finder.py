"""
root_finder.py

to find the roots of mulitvariate function, being aerodynamic torques, which is
a function of rotation speed, pitch and cone angle.

1. load data

"""

# dependencies. SciPy has the ND interpolation function
import numpy as np
from scipy.interpolate import LinearNDInterpolator, griddata
from scipy.optimize import root
import matplotlib.pyplot as plt
import pyvista as pv


def load_data(filename):
    """loads csv data to interpolate.\n
    should be of the form:\n
    Case, Cone Angle, Omega, Pitch Angle, forces (x3), moments (x3).\n
    RETURNS: cone, omega, pitch, mx, my, mz"""

    # csv, first row names, first column is non-numeric so ignore
    data = np.loadtxt(filename, delimiter=',', skiprows=1, usecols=(1,2,3,4,5,
                                                                    6,7,8,9))

    # run variables
    cone = data[:,0]; omega = data[:,1]; pitch = data[:,2]
    variables = data[:,0:3]

    # forces (ignore for now)
    fx = data[:,3]; fy = data[:,4]; fz = data[:,5]

    # moments
    mx = data[:,6]; my = data[:,7]; mz = data[:,8]
    moments = data[:,6:]

    # return cone, omega, pitch, mx, my, mz
    return variables, moments


def plot_moments(filename):
    """to plot moments"""
    variables, moments = load_data(filename)

    cone = variables[:,0]
    omega = variables[:,1]
    pitch = variables[:,2]

    mx = moments[:,0]
    my = moments[:,1]
    mz = moments[:,2]

    ####################################################################

    # x = cone, y = omega, z = pitch
    xmin = min(variables[:,0])
    xmax = max(variables[:,0])

    ymin = min(variables[:,1])
    ymax = max(variables[:,1])

    zmin = min(variables[:,2])
    zmax = max(variables[:,2])


    # setting up grid
    N = 20

    xs = np.linspace(xmin, xmax, N)
    ys = np.linspace(ymin, ymax, N)
    zs = np.linspace(zmin, zmax, N)

    X, Y, Z = np.meshgrid(xs, ys, zs, indexing='ij')

    R = griddata(
            points=(cone, omega, pitch),
            values=mx,
            xi=(X, Y, Z),
            method='linear'
            )            

    R = np.nan_to_num(R, nan=0.0)  # replace nans with zero for safe plotting

    grid = pv.StructuredGrid()
    grid.points = np.column_stack((X.flatten(), Y.flatten(), Z.flatten()))
    grid.dimensions = X.shape
    grid["Residual"] = R.flatten()

    # ------------------------------
    # Plot in full 3D volume
    # ------------------------------
    plotter = pv.Plotter()
    plotter.add_volume(grid, opacity='sigmoid', cmap='viridis')
    plotter.add_axes()
    plotter.show(title="3D Residual Field (‖F(x,y,z)‖)")



    
    ####################################################################



def interpolate(filename):
    """input: filename. perform NDinterpolation"""

    variables, moments = load_data(filename)

    idx = np.argmin(np.linalg.norm(moments, axis=1))
    x0 = variables[idx]

    # --- Step 2: interpolate ---
    Fx = LinearNDInterpolator(variables, moments[:,0])
    Fy = LinearNDInterpolator(variables, moments[:,1])
    Fz = LinearNDInterpolator(variables, moments[:,2])

    def F(v):
        return np.array([Fx(*v), Fy(*v), Fz(*v)])

    # --- Step 3: solve root ---
    sol = root(F, x0)
    print("Root:", sol.x)
    print("Residual:", sol.fun)
    print("Success:", sol.success)
    print("Message:", sol.message)

    print(Fx(*sol.x), Fy(*sol.x), Fz(*sol.x))
    print("Min residual per component:", moments.min(axis=0))
    print("Max residual per component:", moments.max(axis=0))


filename = 'kateForces.csv'
interpolate(filename)
plot_moments(filename)
