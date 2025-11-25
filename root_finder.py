"""
root_finder.py

to find the roots of mulitvariate function, being aerodynamic torques, which is
a function of rotation speed, pitch and cone angle.

1. load data
2. Interpolate data to a continuum
3. root find

"""

# dependencies. SciPy has the ND interpolation function
import numpy as np
from scipy.interpolate import RBFInterpolator, griddata
from scipy.optimize import root, least_squares
import matplotlib.pyplot as plt
from matplotlib.ticker import LinearLocator
import argparse


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
    # setting up subplots
    fig = plt.figure(figsize=(12, 4))

    # text object to tell the time on the video
    time_text = fig.text(0.05, 0.8, "", fontsize=12)

    # ax = fig.add_subplot(111, projection='3d')
    
    ax1 = fig.add_subplot(131, projection='3d')
    ax2 = fig.add_subplot(132, projection='3d')
    ax3 = fig.add_subplot(133, projection='3d')

    # x = cone, y = omega, z = pitch
    xmin = min(variables[:,0])
    xmax = max(variables[:,0])

    ymin = min(variables[:,1])
    ymax = max(variables[:,1])

    zmin = min(variables[:,2])
    zmax = max(variables[:,2])

    # for ax in (ax1, ax2, ax3):
    #     ax.set_xlim(xmin,xmax)
    #     ax.set_ylim(ymin,ymax)
    #     ax.set_zlim(zmin,zmax)

    # ax1.set_title("All over view")
    # ax2.set_title("Top (Z)")
    # ax3.set_title("Side (Y)")
    ####################################################################

    # setting up grid
    N = 20

    xs = np.linspace(xmin, xmax, N)
    ys = np.linspace(ymin, ymax, N)
    zs = np.linspace(zmin, zmax, N)

    X, Y = np.meshgrid(xs, ys)
    Z = zs

    surface = ax.plot_surface(X, Y, Z, cmap='viridis', linewidth=0)

    ax.set_zlim(zmin, zmax)
    ax.zaxis.set_major_locator(LinearLocator(10))
    # A StrMethodFormatter is used automatically
    ax.zaxis.set_major_formatter('{x:.02f}')

    # Add a color bar which maps values to colors.
    fig.colorbar(surface, shrink=0.5, aspect=5)

    plt.show()

    ####################################################################

def interpolate(variables, moments):

    mx = moments[:,0]
    my = moments[:,1]
    mz = moments[:,2]

    mx_rbf = RBFInterpolator(variables, mx,kernel="thin_plate_spline")
    my_rbf = RBFInterpolator(variables, my,kernel="thin_plate_spline")
    mz_rbf = RBFInterpolator(variables, mz,kernel="thin_plate_spline")

    return mx_rbf, my_rbf, mz_rbf



def root_find(filename, moment_target=[0.0, 0.0, 0.0], weights=[1.0,1.0,1.0]):
    """input: filename and target moment\n
    these are the equal and opposite moments from forces due to gravity!\n
    performs radial basis function interpolation to get variables as a 
    continuum and then 
    least squares or
    root()"""

    # loading data into arrays
    moment_target = np.array(moment_target)
 
    variables, moments = load_data(filename)

    # optional argument, (to try and minimise torque about x and y more for example)
    weights = np.array(weights)

    # builds continuous interpolated functions out of look up table
    mx_rbf, my_rbf, mz_rbf = interpolate(variables, moments)

    ###########################################################################
    # functions to call in scipy solvers

    def M_interpolated(x):
        """M(x) evaluated at x, x is column vector of variables"""
        # scipy interpolators require shape 1,3, and return arrays
        v = np.asarray(x).reshape(1, -1)
        # convert array to single scalar value
        return np.array([mx_rbf(v)[0], my_rbf(v)[0], mz_rbf(v)[0]])

    def residual_moment(x):
        return weights * (M_interpolated(x) - moment_target)
    
    ###########################################################################

    # initial guess - the norm value that minimises all torques
    sample_diff_norms = np.linalg.norm((moments - moment_target) * weights, axis=1)
    i0 = int(np.argmin(sample_diff_norms))
    x0 = variables[i0]   

    bounds_min = variables.min(axis=0)
    bounds_max = variables.max(axis=0)

    # need leastsqrs here because the data contains no roots!!
    res_ls = least_squares(residual_moment, x0, bounds=(bounds_min, bounds_max), xtol=1e-12, ftol=1e-12)

    print("#####    Least-squares result:    #####")
    print(" - minimal variables:", res_ls.x)
    print(" - residuals:", res_ls.fun)
    print(" - success:", res_ls.success)


# needs tidying
# fpath = 'OF_data/kateForces.csv'
# fpath2 = 'OF_data/forces_results.csv'
# root_find(argv[0])



if __name__ == '__main__':
    parser = argparse.ArgumentParser(description='''
                                     running root finder code on openfoam data
                                     ''')
    parser.add_argument('filepath', type=str, help='path to data file')

    args = parser.parse_args()

    fpath = args.filepath

    root_find(fpath)
