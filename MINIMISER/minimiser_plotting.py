"""
root_finder.py

to find the roots of mulitvariate function, being aerodynamic torques, 
which is a function of rotation speed, pitch and cone angle.

1. load data
2. Interpolate data to a continuum
3. root find

"""

# dependencies. SciPy has the NDim interpolation function
import numpy as np
from scipy.interpolate import RBFInterpolator, RegularGridInterpolator, griddata, LinearNDInterpolator
from scipy.optimize import root, least_squares
from scipy.spatial.transform import Rotation as R
import matplotlib.pyplot as plt
import argparse

import minimiser as mn


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


def root_find_plot(filename, inertia, cone, pitch, vz, mass=0.0281, g=9.81, epsilon=1.0, cover=1.0, window=1.0):
    """routine for plotting net moments and forces for slice of data
    input: filename and target moment\n
    output: prints minimised values to terminal\n
    these are the equal and opposite moments from forces due to gravity!\n
    performs radial basis function interpolation to get variables as a 
    continuum and then least squares or root()"""

    print("""########################################################################
force scaling factor: """, epsilon)

    # load data
    variables, forces, moments = load_data(filename)

    # for taking random elements of the whole dataset
    rng = np.random.default_rng(seed=42)
    num_cases = variables.shape[0]
    num_samples = int(num_cases * cover)
    print(f"Number of original cases: {num_cases}, number to sample: {num_samples}.")
    indices = rng.choice(num_cases, size=num_samples, replace=False)

    variables = variables[indices,:]
    moments = moments[indices,:]
    fzs = forces[indices, 2]

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

        I2 = rotation_matrix(180, 0, I)

        # I = R @ I0 @ R.T
        I_rot = rotation_matrix(cone, pitch, I2)

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

    fzNet = []
    MxNet = []
    MyNet = []
    MzNet = []
    vzs = []
    omegas = []

    mx_max = max(moments[:,0])
    my_max = max(moments[:,1])
    mz_max = max(moments[:,2])
    fz_max = max(fzs)


    norm_factor = np.array([mx_max, my_max, mz_max, fz_max])

    for i, m in enumerate(moments):

        cfd_res = np.array([m[0], m[1], m[2], epsilon*fzs[i]])

        inertial_res = mass_moment(variables[i,:])

        difference[i] = np.linalg.norm(cfd_res + inertial_res)

        x = variables[i,:]
        if x[0] == cone and x[2] == pitch and x[3] == vz:
            omegas.append(x[1])

            mnet = m + inertial_res[:3]
            MxNet.append(mnet[0])
            MyNet.append(mnet[1])
            MzNet.append(mnet[2])
            fzNet.append(fzs[i]*epsilon+inertial_res[3])

    # index of that guess
    i0 = int(np.argmin(difference))
    x0 = variables[i0]

    min_diffs = np.argsort(difference)[:4]
    print("min diff", variables[min_diffs[1]],
variables[min_diffs[2]],
variables[min_diffs[3]])
    
    fzs = np.array(fzs)
    omegas = np.array(omegas)
    vzs = np.array(vzs)

    MxNet = np.array(MxNet)    # Net moments in x
    MyNet = np.array(MyNet)    # Net moments in y
    MzNet = np.array(MzNet)    # Net moments in z
    fzNet = np.array(fzNet)    # Net forces

    sorted_idx = np.argsort(omegas)
    sorted_omegas = omegas[sorted_idx]

    plt.figure(figsize=(7,5))
    plt.title("Net Force and Moments", fontsize=17)
    plt.ylabel("Residual Force and Moments", fontsize=15)
    plt.xlabel("Rotation Rate (rpm)", fontsize=15)
    plt.tick_params(labelsize=13)
    plt.grid()
    
    plt.plot(sorted_omegas, MxNet[sorted_idx], label=f'$M_x$', marker='o')
    
    plt.plot(sorted_omegas, MyNet[sorted_idx], label='$M_y$', marker='s')
    plt.plot(sorted_omegas, MzNet[sorted_idx], label='$M_z$', marker='v')
    # plt.plot(omegas, Mnet, label=f'$M_x$', marker='o')

    plt.plot(sorted_omegas, fzNet[sorted_idx],  label='$F_z$', marker='x', color='purple')

    plt.axhline(y=0, ls=':', color='red')

    plt.legend(fontsize=13)

    plt.text(0.7, 0.1,f"Cone=5$^{{\circ}}$\nPitch=-4$^{{\circ}}$\n$v_z$=3.25 m/s", transform=plt.gca().transAxes, fontsize=12, bbox={'facecolor': 'white', 'alpha': 0.6})

    plt.savefig("../figures/netZeroRoot.png", bbox_inches='tight', dpi=300)
    plt.close()

###############################################################################
    # for 3D plot

    # from matplotlib import cm
    # from matplotlib.colors import Normalize

    # # Calculate omegas/vzs
    # x_vals = omegas / vzs  # x-axis: omegas/vzs
    # y_vals = pitches  # y-axis: pitch angles
    # z_vals = np.array(sumDiff)  # z-axis: sumDiff
    # vz_vals = np.unique(vzs)  # Unique vz values

    # # Normalize vz values for coloring
    # norm = Normalize(vmin=vz_vals.min(), vmax=vz_vals.max())
    # cmap = cm.viridis  # Choose a colormap

    # # Create the 3D plot
    # fig = plt.figure(figsize=(7, 5))
    # ax = fig.add_subplot(111, projection='3d')
    # fig.tight_layout()

    # # Loop through each unique vz value
    # for vz in vz_vals:
    #     # Filter data for the current vz
    #     mask = vzs == vz
    #     x_filtered = x_vals[mask]
    #     y_filtered = y_vals[mask]
    #     z_filtered = z_vals[mask]

    #     # Create a grid for omegas/vzs and pitch angles
    #     x_grid = np.linspace(x_filtered.min(), x_filtered.max(), 100)
    #     y_grid = np.linspace(y_filtered.min(), y_filtered.max(), 100)
    #     x_mesh, y_mesh = np.meshgrid(x_grid, y_grid)

    #     # Interpolate sumDiff values for the grid
    #     interpolator = LinearNDInterpolator(list(zip(x_filtered, y_filtered)), z_filtered)
    #     z_mesh = interpolator(x_mesh, y_mesh)

    #     # Get the color for the current vz
    #     color = cmap(norm(vz))

    #     # Create a 2D array of colors matching the shape of z_mesh
    #     facecolors = np.full((*z_mesh.shape,4), color)

    #     # Plot the surface for the current vz
    #     surf = ax.plot_surface(x_mesh, y_mesh, z_mesh, facecolors=facecolors, edgecolor='none', alpha=0.8)

    # # # Add a zero-level surface
    # # x_vals = np.linspace(ax.get_xlim()[0], ax.get_xlim()[1], 100)
    # # y_vals = np.linspace(ax.get_ylim()[0], ax.get_ylim()[1], 100)
    # # X, Y = np.meshgrid(x_vals, y_vals)
    # # Z = np.zeros_like(X)  # Zero-level surface

    # # # Plot the zero-level surface
    # # ax.plot_surface(X, Y, Z, alpha=0.3, color='red', label='Net Zero')

    # # Add a color bar to show the vz scale
    # mappable = cm.ScalarMappable(norm=norm, cmap=cmap)
    # mappable.set_array(vz_vals)
    # cbar = plt.colorbar(mappable, ax=ax, shrink=0.5, aspect=10)
    # cbar.set_label('$v_z$ (m/s)', fontsize=15)
    # cbar.ax.tick_params(labelsize=13) 

    # # Set axis labels
    # ax.set_xlabel('Rotation Speed / Descent Velocity (rpm / (m/s))', fontsize=15, labelpad=10)
    # ax.set_ylabel('Pitch Angle (degrees)', fontsize=15, labelpad=10)
    # ax.set_zlabel('Residual From Net Zero', fontsize=15, labelpad=10)
    # ax.tick_params(labelsize=13)

    # # Set title
    # ax.set_title('Net Combined Moment and Force', fontsize=17, pad=-20)
###############################################################################

    # plt.show()
    # x0 = [4, 420, -1, 3]  #### can set the first guess here

    mm = mass_moment(x0)

    print("""variables: cone rpm pitch velocity
init guess""", x0)
    print(f"mass moment about first guess: {mm}")
    print(f"cfd moment about first guess:  {moments[i0]} {fzs[i0]}")
    print(f"Net moments: {mm[:3]+moments[i0]}\nNet force: {fzs[i0]-mass*g}")

    # this is array type - min and max angles and omegas
    bounds_min = variables.min(axis=0)
    bounds_max = variables.max(axis=0)

    # bounds_min = [ -10.,   0.,   -20.,    0.5]
    # bounds_max = [ 30.,  800.,   15.,    7.5]

    print("min", bounds_min)
    print("max", bounds_max)

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
    print(" - residuals:", result_ls.fun)
    print(" - success:", result_ls.success)
    
    return min_variables


def plot_continuum(filename, inertia, mass, g, epsilon, cover=1, plot=True):
    # fixed parameters to plot 2d graph
    vz = 3
    cone_angle = 3
    pitch_angle = -1
    nPoints = 100
    # load data
    variables, forces, moments = load_data(filename)
    # and apply windowing about the Fz = mg cases
    # variables, forces, moments = windowing(variables, forces, moments, mass, g)
    orig_vars = variables
    orig_forces = forces
    orig_moments = moments

    # for taking random elements of the whole dataset
    # rng = np.random.default_rng(seed=42)
    # num_cases = variables.shape[0]
    # num_samples = int(num_cases * cover)
    # indices = rng.choice(num_cases, size=num_samples, replace=False)

    # variables = variables[indices,:]
    # moments = moments[indices,:]
    # fzs = forces[indices, 2]

    variables = np.array([x for i, x in enumerate(orig_vars) if (i+1) % cover != 0])
    moments = np.array([x for i, x in enumerate(orig_moments) if (i+1) % cover != 0])
    fzs = np.array([x[2] for i, x in enumerate(orig_forces) if (i+1) % cover != 0])

    # builds continuous interpolated functions out of look up table
    mx_rbf, my_rbf, mz_rbf, fz_rbf = interpolate(variables, moments, fzs*epsilon)

    # loading inertia about principle axes into 3x3 array
    I = np.zeros(shape=(3,3))

    I[0,0] = inertia[0]
    I[1,1] = inertia[1]
    I[2,2] = inertia[2]

    def cfd_moment(x):
        """M(x) evaluated at x, x is column vector of variables"""
        moment = np.array([mx_rbf([x])[0], my_rbf([x])[0], mz_rbf([x])[0], fz_rbf([x])[0]])
        return moment

    def mass_moment(x):
        """input: variables [cone, rpm, pitch, vz]\n
        output: the target moments and force to be equal and opposite to
        [moment[0], moment[1], moment[2], weight]"""
        cone = x[0]
        pitch = x[2]
        omega = np.array([0, 0, x[1]*np.pi/30])

        # I = R @ I0 @ R.T
        I_rot = rotation_matrix(cone, pitch, I)

        # L = I @ omega
        L = np.matmul(I_rot, omega)

        # M = omega x L
        moment = np.cross(omega, L)

        weight = -mass * g * epsilon

        result = np.array([moment[0], moment[1], moment[2], weight])
        return result

    N = len(variables[:,1])

    cData = variables[:,0]
    oData = variables[:,1]
    pData = variables[:,2]
    vData = variables[:,3]

    mass_moments = []
    omegas = []
    table_moments = []
    fz = []

    # searching for matching cases in the skipped and interpolated data
    for i in range(N):
        # cases.append([cones[i], omegas[i], pitches[i], vzs[i]])
        # try for changing only omega
        x = variables[i,:]
        if x[0] == cone_angle and x[2] == pitch_angle and x[3] == vz:
            omegas.append(oData[i])
            mass_moments.append(mass_moment([cone_angle, oData[i], pitch_angle, vz]))
            table_moments.append(moments[i,:])
            fz.append(epsilon*forces[i,2])

    # for the whole, unskipped data set comparison
    orig_omegas = []
    orig_moments2 = []
    orig_fz = []

    for i, x in enumerate(orig_vars):
        if x[0] == cone_angle and x[2] == pitch_angle and x[3] == vz:
            orig_omegas.append(x[1])
            orig_moments2.append(orig_moments[i,:])
            orig_fz.append(epsilon*orig_forces[i,2])

    orig_moments2 = np.array(orig_moments2)

    # cfd_moments = [cfd_moment(x) for x in cases]
    # cfd_moments = np.array(cfd_moments)
    mass_moments = np.array(mass_moments)
    omegas = np.array(omegas)
    table_moments = np.array(table_moments)
    fz = np.array(fz)

    omega_interp = np.linspace(np.min(omegas), np.max(omegas), nPoints)

    cfd_moments = [cfd_moment([cone_angle, x, pitch_angle, vz]) for x in omega_interp]
    cfd_moments = np.array(cfd_moments)

    # M_mass = [-0.00158594,  0.01497816,  0]

    # plotting
    if plot:
        coords = ['X-Moment', 'Y-Moment', 'Z-Moment', 'Drag Force']

        for i in range(3):
            fig, axs = plt.subplots(figsize=(7,5))

            axs.set_title(f"{coords[i]} Interpolated, cone = {cone_angle}, pitch = {pitch_angle}, vz = {vz}", fontsize=17)

            axs.set_ylabel("Moment (Nm)", fontsize=15)
            axs.set_xlabel("Rotation Speed (rpm)", fontsize=15)

            axs.plot(omega_interp, cfd_moments[:,i], label='Interpolation')

            axs.scatter(orig_omegas, orig_moments2[:,i], marker='x', color='purple', label='CFD results (all)')
            axs.scatter(omegas, table_moments[:,i], color='r', label='CFD results (incl.)')
            # axes[i].plot(omegas,-mass_moments[:,i], ls='--')
            axs.tick_params(labelsize=13)
            axs.grid()
            axs.legend(fontsize=13)
            plt.savefig(f"../figures/{coords[i]}_interp.png", bbox_inches='tight', dpi=300)
            plt.close()

        fig, axs = plt.subplots(figsize=(7,5))
        # Force in Z direction
        fig.suptitle(f"Drag Force, cone = {cone_angle}, pitch = {pitch_angle}, vz = {vz}", fontsize=16)
        axs.plot(omega_interp, cfd_moments[:,3])
        axs.scatter(orig_omegas, orig_fz, color="r")

        axs.set_xlabel("Rotation Speed (rpm)", fontsize=14)
        axs.set_ylabel(f"Drag Force ({epsilon}N)", fontsize=14)
        axs.grid()
        plt.savefig(f"../figures/forceZ_interp.png", bbox_inches='tight', dpi=300)
        plt.close()

    ###########################################################################
    # error analysis

    # interpolated moments for each true data point
    # BUT the interpolation is still built on the sparse data (i think)
    interp_moments = [cfd_moment([cone_angle, x, pitch_angle, vz]) for x in orig_omegas]
    interp_moments = np.array(interp_moments)

    residual_x = np.abs( (orig_moments2[:,0] - interp_moments[:,0]) / orig_moments2[:,0] * 100 )
    residual_y = np.abs( (orig_moments2[:,1] - interp_moments[:,1]) / orig_moments2[:,1] * 100 )
    residual_z = np.abs( (orig_moments2[:,2] - interp_moments[:,2]) / orig_moments2[:,2] * 100 )
    return residual_x, residual_y, residual_z


def test_interp_thinning(filename, epsilon, inertia, mass, g, N):

    residuals = []
    mean_res = []
    results = []
    unsampled = []

    scale = 1.0

    for i in range(N):
        # orig_omegas gets written over but thats okay it should always be the same
        res_x, res_y, res_z = plot_continuum(filename, inertia, mass, g, epsilon, cover=(N-i*scale), plot=False)
        # residuals.append(residual)
        print("res_x", len(res_x))

        residual_x = [x for x in res_x if str(x) != 'nan']
        residual_y = [x for x in res_y if str(x) != 'nan']
        residual_z = [x for x in res_z if str(x) != 'nan']
        mean_res.append((np.mean(residual_x), np.mean(residual_y), np.mean(residual_z)))
        undef = len(res_x) - len(residual_x)
        unsampled.append(undef)

        result = mn.root_find(filename, inertia, mass=mass, epsilon=epsilon, cover=(N-i*scale))
        results.append(result)

    print(mean_res)

    results = np.array(results)
    mean_res = np.array(mean_res)
    # mean_res = [(np.float64(22.912606090744628), np.float64(1.6731643524406663), np.float64(82.13330253799352)), (np.float64(24.382318594104355), np.float64(1.6958424226456974), np.float64(63.13363345208543)), (np.float64(24.382318594104326), np.float64(1.6958424226457134), np.float64(63.133633452085384)), (np.float64(40.37349708407162), np.float64(2.9107038329001678), np.float64(116.39960779056852)), (np.float64(28.989882798860272), np.float64(3.182364559457358), np.float64(78.18733823587041)), (np.float64(41.92963899841839), np.float64(5.015154047072699), np.float64(92.89627111560496)), (np.float64(34.3792312943962), np.float64(3.9224076756855486), np.float64(74.9940885281133)), (np.float64(46.4783778235721), np.float64(6.405534005730114), np.float64(105.9579542052004)), (np.float64(36.13920666872251), np.float64(6.564666760053472), np.float64(96.3493709977765)), (np.float64(76.91144466846316), np.float64(10.842487806366002), np.float64(155.7233513881097))]

    plt.figure(figsize=(7,5))
    plt.title("Actual vs Interpolated Moments After Data Removal", fontsize=17)
    plt.ylabel("Difference in Mean Values (%)", fontsize=15)
    plt.xlabel("Percentage Data Removed (%)", fontsize=15)
    plt.tick_params(labelsize=13)
    plt.grid()
    plt.plot([(i*scale)/N*100 for i in range(N)], [x[0] for x in mean_res], marker='o', label="Mx")
    plt.plot([(i*scale)/N*100 for i in range(N)], [x[1] for x in mean_res], marker='s', label="My")
    plt.plot([(i*scale)/N*100 for i in range(N)], [x[2] for x in mean_res], marker='v', label="Mz")

    plt.legend(fontsize=12)

    plt.savefig("../figures/dataRemovalMY.png", bbox_inches='tight', dpi=300)
    plt.close()


    fig, ax1 = plt.subplots(figsize=(9,6))

    ax1.set_title("Actual vs Interpolated Moments About x", fontsize=17)

    color1 = "steelblue"
    ax1.set_ylabel("Difference in Mean Values (%)", fontsize=15)
    ax1.set_xlabel("Data Removed (%)", fontsize=15)
    ax1.plot([(i*scale)/N*100 for i in range(N)], [x[0] for x in mean_res], marker='o', label="Mean")
    ax1.tick_params(axis='y', labelcolor=color1, labelsize=13)
    ax1.set_xscale("log")

    ax1.grid()
    ax2 = ax1.twinx()

    color2 = "orange"
    ax2.set_ylabel("Number of Original Points Not Covered By Interpolation", color=color2, fontsize=15)
    ax2.plot([(i*scale)/N*100 for i in range(N)], unsampled, color=color2, marker='s', label="No. Unsampled")
    ax2.tick_params(axis='y', labelcolor=color2, labelsize=13)

    ax1.legend(fontsize=13)

    fig.tight_layout()
    plt.savefig("../figures/dataRemoval-twinx.png", bbox_inches='tight', dpi=300)
    plt.close()

    print(results)

    init_result = results[0,:]

    print(init_result)
    # results = results[1:,:]

    plt.figure(figsize=(7,5))
    plt.title("Steady-State Solution After Data Removal", fontsize=17)
    plt.ylabel("Variable with Data Removed / Without", fontsize=15)
    plt.xlabel("Percentage Data Removed (%)", fontsize=15)
    plt.tick_params(labelsize=13)

    plt.plot([(i*scale)/N*100 for i in range(N)], [x[0]/init_result[0] for x in results], marker='o', label='Cone')
    plt.plot([(i*scale)/N*100 for i in range(N)], [x[1]/init_result[1] for x in results], marker='x', label='RPM')
    plt.plot([(i*scale)/N*100 for i in range(N)], [x[2]/init_result[2] for x in results], marker='s', label='Pitch')
    plt.plot([(i*scale)/N*100 for i in range(N)], [x[3]/init_result[3] for x in results], marker='v', label='$v_z$')

    plt.grid()
    plt.legend(fontsize=13)
    plt.savefig("../figures/dataRemovalSS.png", bbox_inches='tight', dpi=300)
    plt.close()

    # residuals = np.array(residuals)

    # for i in range(N):
    #     plt.plot(orig_omegas, residuals[i,:])
    # plt.show()


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

    dst = ''  # where to save figures to

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

    mn.root_find(fpath, inertia=inertia, mass=mass, g=g, epsilon=epsilon, 
             cover=2, window=window)

    plot_continuum(fpath, inertia, mass, g, epsilon, cover=4)

    test_interp_thinning(fpath, epsilon, inertia, mass, g, 10)

    # test_thinning(fpath)


