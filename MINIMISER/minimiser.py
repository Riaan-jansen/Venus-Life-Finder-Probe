"""
minimiser.py

to find the minima of mulitvariate function, being total moments, and
lift which are a function of rotation speed, pitch and cone angle.

1. load data (optionally windowing around F=mg and scaling force by a factor)
2. Interpolate data to a continuum
3. find the minimum of the combined moments and force from the data table
4. find the minimum from the interpolated functions

"""

# dependencies. SciPy has the NDim interpolation function
import numpy as np
from scipy.interpolate import LinearNDInterpolator
from scipy.optimize import least_squares
from scipy.spatial.transform import Rotation as R
import matplotlib.pyplot as plt
import argparse

import minimiser_plotting as mpl


def load_data(filename):
    """loads csv data to interpolate.\n
    should be of the form:\n
    Case, Cone Angle, Omega, Pitch Angle, vz, forces (x3), moments (x3).\n
    Returns: variables, forces, moments
    Requires specific ordering of csv file."""

    # csv, first row names, first column is non-numeric so ignore
    data = np.loadtxt(filename, delimiter=',', skiprows=1, 
                      usecols=(1,2,3,4,5,6,7,8,9,10), comments='#')

    # extract run variables (angles, rpm, vz) and forces and moments
    variables = data[:,0:4]
    forces = data[:,4:7]
    moments = data[:,7:]

    return variables, forces, moments


def windowing(variables, forces, moments, mass, g, window):
    """to search only around Fz = mg cases\n
    input: variables, forces, moments, mass, g\n
    output: new variables, forces and moments"""

    new_vars = []
    new_forces = []
    new_moments = []

    # actual probe weight
    weight = mass*g

    for i, f in enumerate(forces):
        fz = f

        # chop the data to only look at near terminal velocity runs
        if np.abs(fz-weight)/weight <= window:
            new_vars.append([x for x in variables[i,:]])
            # new_forces.append([x for x in forces[i]])
            new_forces.append(f)
            new_moments.append([x for x in moments[i,:]])

    new_vars = np.array(new_vars)
    new_forces = np.array(new_forces)
    new_moments = np.array(new_moments)

    return new_vars, new_forces, new_moments


def interpolate(variables, moments, fz):
    """linear interpolator"""
    fill = 1e6  # penalty for going outside the data bounds
    mx = LinearNDInterpolator(variables, moments[:, 0], fill_value=fill)
    my = LinearNDInterpolator(variables, moments[:, 1], fill_value=fill)
    mz = LinearNDInterpolator(variables, moments[:, 2], fill_value=fill)
    fz = LinearNDInterpolator(variables, fz,            fill_value=fill)

    return mx, my, mz, fz


def rotation_matrix(cone, pitch, I):
    """ rotating inertia tensor by cone and then pitch (reverse of how 
    we arrived there.
    I_new = R.I.R^T 
    credit to James."""

    rot = R.from_euler(seq='xy', angles=[pitch, -cone], degrees=True)

    rot_mat = rot.as_matrix()
    t_rot_mat = np.transpose(rot_mat)

    I_rot = np.matmul(I, rot_mat)
    I_final = np.matmul(t_rot_mat, I_rot)

    return I_final


def matmul(a, b):
    """input: a (3x3), b (3)\n
    some reason standard methods wont give a 3x3 matrix"""
    result = np.zeros((3,3))

    # iterate over rows
    for i in range(3):
        result[0,i] = (a[0,i] + a[1,i] + a[2,i]) * b[0]
        result[1,i] = (a[0,i] + a[1,i] + a[2,i]) * b[1]
        result[2,i] = (a[0,i] + a[1,i] + a[2,i]) * b[2]

    return result


def root_find(filename, inertia, mass=0.0281, g=9.81, epsilon=1.0, cover=1, window=1.0):
    """input: filename and target moment\n
    output: prints minimised values to terminal\n
    these are the equal and opposite moments from forces due to gravity!\n
    performs radial basis function interpolation to get variables as a 
    continuum and then least squares or root()"""

    print("""########################################################################
force scaling factor: """, epsilon)

    # load data
    variables, forces, moments = load_data(filename)

    # for taking random elements of the whole dataset
    # rng = np.random.default_rng(seed=42)
    # num_cases = variables.shape[0]
    # num_samples = int(num_cases * cover)
    # print(f"Number of original cases: {num_cases}, number to sample: {num_samples}.")
    # indices = rng.choice(num_cases, size=num_samples, replace=False)

    # variables = variables[indices,:]
    # moments = moments[indices,:]
    # fzs = forces[indices, 2]
    orig_vars = variables
    orig_forces = forces
    orig_moments = moments

    print(np.shape(orig_moments))

    variables = np.array([x for i, x in enumerate(orig_vars) if (i+1) % cover != 0])
    moments = np.array([x for i, x in enumerate(orig_moments) if (i+1) % cover != 0])
    fzs = np.array([x[2] for i, x in enumerate(orig_forces) if (i+1) % cover != 0])

    # variables = np.array(variables)
    # moments = np.array(moments)
    # fzs = np.array(fzs)
    print(np.shape(moments))


    # print("number of cases: ", len(variables))
    # # and apply windowing about the Fz = mg cases
    variables, fzs, moments = windowing(variables, fzs, moments, mass, g, window)

    # # builds continuous interpolated functions out of look up table
    mx_rbf, my_rbf, mz_rbf, fz_rbf = interpolate(variables, moments, fzs*epsilon)

    print("no. cases after window: ", variables.shape[0])

    # loading inertia about principle axes into 3x3 array
    I = np.zeros(shape=(3,3))

    I[0,0] = inertia[0]
    I[1,1] = inertia[1]
    I[2,2] = inertia[2]

###########################################################################
    # functions to call in scipy solvers

    def cfd_moment(x):
        """M(x) evaluated at x, x is column vector of variables"""
        moment = np.array([mx_rbf([x])[0], my_rbf([x])[0], mz_rbf([x])[0], fz_rbf([x])[0]])
        return moment

    def net_moment(x):
        """( aerodynamic moment - moment due to mass ) * weighting"""
        # net_M = np.linalg.norm(cfd_moment(x) + mass_moment(x))
        net_M = (cfd_moment(x) + mass_moment(x))

        return net_M

    def mass_moment(x):
        """input: variables [cone, rpm, pitch, vz]\n
        output: the target moments and force to be equal and opposite to
        [moment[0], moment[1], moment[2], weight]"""
        cone = x[0]
        pitch = x[2]
        omega = np.array([0, 0, x[1]*np.pi/30])
        # vz is essentially a dummy variable here. has no effect on target Fz
        # but the answer still returns the best guess as a result of varying vz
        # over different runs and the combination of pitch and cone etc that
        # lead to minimised moments
        # vz = x[3]

        # I2 = rotation_matrix(180, 0, I)

        # I = R @ I0 @ R.T
        I_rot = rotation_matrix(cone, pitch, I)

        # L = I @ omega
        L = np.matmul(I_rot, omega)

        # M = omega x L
        moment = np.cross(omega, L)

        weight = -mass * g * epsilon

        result = np.array([moment[0], moment[1], moment[2], weight])
        return result

    ###########################################################################

    # initial guess - the norm value that minimises all torques
    difference = np.zeros(len(moments[:,0]))

    # this might be important, technically cant add forces and moments without
    # some factor to account for the different dimensions
    mx_max = max(moments[:,0])
    my_max = max(moments[:,1])
    mz_max = max(moments[:,2])
    fz_max = max(fzs)
    norm_factor = np.array([mx_max, my_max, mz_max, fz_max])

    for i, m in enumerate(moments):

        cfd_res = np.array([m[0], m[1], m[2], epsilon*fzs[i]])

        inertial_res = mass_moment(variables[i,:])

        difference[i] = np.linalg.norm(cfd_res + inertial_res)

    # index of that guess
    i0 = int(np.argmin(difference))
    x0 = variables[i0]

#     min_diffs = np.argsort(difference)[:4]
#     print("Next best guesses: ", variables[min_diffs[1]],
# variables[min_diffs[2]],
# variables[min_diffs[3]])

    # x0 = [4, 420, -1, 3]  #### can set the first guess here

    mm = mass_moment(x0)

    print(f"""variables: cone rpm pitch velocity
initial guess: {x0}
mass moment about first guess: {mm}
cfd moment about first guess:  {moments[i0]} {fzs[i0]}
Net moments: {mm[:3]+moments[i0]}
Net force: {fzs[i0]-mass*g}
""")

    # this is array type - min and max angles and omegas
    bounds_min = variables.min(axis=0)
    bounds_max = variables.max(axis=0)

    # bounds_min = [ -10.,   0.,   -20.,    0.5]
    # bounds_max = [ 30.,  800.,   15.,    7.5]

    print("Data min: ", bounds_min)
    print("Data max: ", bounds_max)

    # need leastsqrs here because the data contains no roots!!
    result_ls = least_squares(net_moment, x0, bounds=(bounds_min, bounds_max), 
                              max_nfev=1000, method='trf', x_scale=[1,1,1,1])
    min_variables = result_ls.x

    # potentially to then run more simulations around that area of interest.
    print(f"""
######################    Minimisation result:    ######################
 - cone:  {min_variables[0]:.3f}
 - omega: {min_variables[1]:.3f}
 - pitch: {min_variables[2]:.3f}
 - v_t:   {min_variables[3]:.3f}
########################################################################
 - M_mass:{mass_moment(min_variables)}
 - M_cfd :{cfd_moment(min_variables)}
 - M_net :{mass_moment(min_variables)+cfd_moment(min_variables)}
########################################################################""")
    # print(" - residuals:", result_ls.fun)
    print(" - success:", result_ls.success)
    
    return min_variables


if __name__ == '__main__':
    parser = argparse.ArgumentParser(description='''
                                     running root finder code on openfoam data
                                     ''')
    parser.add_argument('filepath', type=str, help='path to data file')

    parser.add_argument('epsilon', type=float, nargs='?', default=1.0, 
                        help='scaling on force')

    parser.add_argument('window', type=float, nargs='?', default=1.0, 
                        help='chose data only around required force by <window> factor')

    args = parser.parse_args()

    fpath = args.filepath
    epsilon = args.epsilon
    window = args.window

    g = 9.8067  # gravitational acceleration. Venus is 8.87 m/s2

    mass = 0.0281         # standard probe
    mass_nut = 0.0383     # added nut mass
    mass_wing = 0.0468    # added nut and wing mass

    inertia=[0.000013713365031784028, 0.0000898853777707176, 
             0.00010289524930712706]
    inertia_nut=[0.0000147227, 0.000103692, 0.000117276]
    inertia_wing=[0.0000208940, 0.000150992, 0.000170734]

    # inertia = [x*1.14 for x in inertia]

    root_find(fpath, inertia=inertia, mass=mass, g=g, epsilon=epsilon, 
              window=window)

# root_find_plot(filename, inertia, cone, pitch, vz, mass=0.0281, g=9.81, epsilon=1.0, cover=1.0, window=1.0)
    mpl.root_find_plot(fpath, inertia, 5, -4, 3.25, mass=mass, g=g, 
                       epsilon=epsilon, window=window)

    # mpl.plot_continuum(fpath, inertia, mass, g, epsilon, cover=0.5)

    # mpl.test_interp_thinning(fpath, epsilon, inertia, mass, g, 7)
