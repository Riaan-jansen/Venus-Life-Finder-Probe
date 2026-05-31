"""
for plotting:
- fz vz scaling
- 
"""

import sys
import numpy as np
import matplotlib.pyplot as plt
from scipy.stats import linregress
from scipy.interpolate import LinearNDInterpolator

from minimiser import windowing
import minimiser as mini


def mass_moment(x, mass, g, I):
    """input: variables [cone, rpm, pitch, vz]\n
    output: the target moments and force to be equal and opposite to
    [moment[0], moment[1], moment[2], weight]"""
    cone = x[0]
    pitch = x[2]
    omega = np.array([0, 0, x[1]*np.pi/30])

    # I = R @ I0 @ R.T
    I_rot = mini.rotation_matrix(cone, pitch, I)

    # L = I @ omega
    L = np.matmul(I_rot, omega)

    # M = omega x L
    moment = np.cross(omega, L)

    weight = -mass * g

    result = np.array([moment[0], moment[1], moment[2], weight])
    return result


def load_data(filename):
    """loads csv data to interpolate. should be of the form:\n
    Case, Cone Angle, Omega, Pitch Angle, vz, forces (x3), moments (x3).\n
    RETURNS: variables, forces, moments"""

    # csv, first row names, first column is non-numeric so ignore
    data = np.loadtxt(filename, delimiter=',', skiprows=1, 
                      usecols=(1,2,3,4,5,6,7,8,9,10), comments='#')

    # extract run variables (angles, rpm, vz) and forces and moments
    variables = data[:,0:4]
    forces = data[:,4:7]
    moments = data[:,7:]

    return variables, forces, moments


def drop_test_result():

    # Data from table
    models = ['Control', 'Nut+', 'Wing+']
    v_up = [2.90, 3.38, 4.83]
    err_up = [0.27, 0.30, 0.44]
    v_down = [2.76, 2.92, 4.15]
    err_down = [0.25, 0.27, 0.38]

    # Create plot
    plt.figure(figsize=(8, 6))

    # Plot "Up" data
    plt.errorbar(models, v_up, yerr=err_up, fmt='-o', capsize=5, label='Up', color='royalblue')

    # Plot "Down" data
    plt.errorbar(models, v_down, yerr=err_down, fmt='-o', capsize=5, label='Down', color='crimson')

    # Formatting
    plt.xlabel('Model')
    plt.ylabel('Terminal Velocity (m/s)')
    plt.title('Terminal Velocity by Model and Direction')
    plt.legend()
    plt.grid(axis='y', linestyle='--', alpha=0.7)

    plt.show()


def plot_fz_vz(filename, cone, pitch, omegas):
    """
    Plot Fz (lift force) vs Vz^2 (descent velocity squared) for different omega values.
    """

    # Load data from the file
    variables, forces, moments = load_data(filename)

    # Initialize lists to store Fz and Vz for each omega
    fz_data = {omega: [] for omega in omegas}
    vz_data = {omega: [] for omega in omegas}

    # Iterate through the data and group by omega
    for i, x in enumerate(variables):
        for omega in omegas:
            if x[0] == cone and x[2] == pitch and x[1] == omega and x[3] >= 2:
                fz_data[omega].append(forces[i, 2])  # Fz (lift force)
                vz_data[omega].append(x[3])         # Vz (descent velocity)

    # Plot Fz vs Vz^2 for each omega
    plt.figure(figsize=(7,5))
    plt.title("Lift Force Scaling with Descent Velocity", fontsize=17)
    plt.xlabel("$v_z^2$ $(m^2s^{-2})$", fontsize=15)
    plt.ylabel("$F_z$ (N)", fontsize=15)
    plt.tick_params(labelsize=13)

    for omega in omegas:
        vz_squared = np.array(vz_data[omega]) ** 2
        fz = np.array(fz_data[omega])
        plt.scatter(vz_squared, fz)  # , label=f"{omega} rpm")
        # Fit a linear regression line
        if len(vz_squared) > 1:  # Ensure there are enough points to fit
            slope, intercept, r_value, _, _ = linregress(vz_squared, fz)
            fz_fit = slope * vz_squared + intercept
            plt.plot(vz_squared, fz_fit, ls="-", label=f"{omega} rpm\n$R^2$={r_value**2:.4f}")


    plt.text(0.7, 0.05,f"Cone={cone}\nPitch={pitch}", transform=plt.gca().transAxes, fontsize=12, bbox={'facecolor': 'white', 'alpha': 0.6})

    plt.grid()
    plt.legend(fontsize=13)
    plt.savefig(dst+"fz-vz.png", bbox_inches="tight", dpi=300)
    plt.close()


def thin_data(variables, moments, fz, n):

    # keeping every nth point
    train_idx = np.arange(0, len(variables), n)
    
    mask = np.ones(len(variables), dtype=bool)
    mask[train_idx] = False
    test_idx = np.where(mask)[0]

    X_train = variables[train_idx]
    X_test  = variables[test_idx]

    fz_train = fz[train_idx]
    fz_test  = fz[test_idx]

    interp = LinearNDInterpolator(X_train, fz_train)

    fz_pred = interp(X_test)

    valid = ~np.isnan(fz_pred)

    coverage = np.mean(valid)

    try:
        rmse = np.sqrt(np.mean((fz_test[valid] - fz_pred[valid])**2))
    except RuntimeWarning:
        rmse = 0

    return coverage, rmse


def test_thinning(filename):

    # load data
    variables, forces, moments = load_data(filename)

    fz = forces[:,2]

    results = []

    # skipN = [1, 2, 3, 4, 5, 6, 7, 8]
    skipN = [2,4,8,16,32,64,128]
    for n in skipN:
        coverage, rmse = thin_data(variables, moments, fz, n)
        results.append((n, coverage, rmse))

    n, cov, err = zip(*results)

    fig, ax1 = plt.subplots(figsize=(9,6))

    fig.suptitle("Interpolation Performance Against Data Reduction", fontsize=15)

    color1 = "steelblue"
    ax1.set_xlabel("n-th data sampled (n)", fontsize=14)
    ax1.set_ylabel("Coverage % (not nan)", color=color1, fontsize=14)
    ax1.plot(n, np.array(cov)*100, marker='o')
    ax1.tick_params(axis='y', labelcolor=color1, labelsize=12)
    ax1.set_xscale("log")

    ax1.grid()
    ax2 = ax1.twinx()

    color2 = "orange"
    ax2.set_ylabel("Error (rms)", color=color2, fontsize=14)
    ax2.plot(n, err, color=color2, marker='s')
    ax2.tick_params(axis='y', labelcolor=color2, labelsize=12)
    ax2.set_ylim(bottom=min(err)*0.75,top=max(err)*1.25)

    fig.tight_layout()
    plt.savefig("../figures/interp-coverage.png", bbox_inches='tight', dpi=300)
    plt.close()


def plot_pitch_Mz(filename, mass, g):

    cone = 5

    variables, forces, moments = load_data(filename)

    # pitch = [x[2] for x in variables]
    # mz = [m[2] for m in moments]
    # fz = [f[2] for f in forces]
    # vz = [x[3] for x in variables]

    plt.figure(figsize=(7,5))
    plt.title("Pitch Angle vs Moment About z", fontsize=17)
    plt.ylabel("Moment (Nm)", fontsize=15)
    plt.xlabel("Pitch Angle (degrees)", fontsize=15)
    plt.tick_params(labelsize=13)
    plt.grid()


    mz = []
    pitches = []
    vzs = []
    omegas = []

    vz = 3

    # Iterate through the data and group by omega
    for i, x in enumerate(variables):

        # if x[0] == cone and x[2] == pitch:
        # if x[2] == pitch:
        if x[0] == cone and x[3] == vz:
            omegas.append(x[1])
            mz.append(moments[i, 2])  # Fz (lift force)
            pitches.append(x[2])         # Vz (descent velocity)
            vzs.append(x[3])

    vzs = np.array(vzs)
    omegas = np.array(omegas)

    scatter = plt.scatter(pitches, mz, c=omegas)
    cbar = plt.colorbar(scatter)
    cbar.set_label( label="Rotation Rate (rpm)", size=15)

    plt.text(0.7, 0.8,f"Cone={cone}$^{{\circ}}$\n$v_z$={vz} m/s", transform=plt.gca().transAxes, fontsize=13, bbox={'facecolor': 'white', 'alpha': 0.6})

    plt.savefig(dst+"pitchVsMz-Vz.png", bbox_inches='tight', dpi=300)
    plt.close()


def plot_omega_Mx(filename, mass, g):

    variables, forces, moments = load_data(filename)

    plt.figure(figsize=(9,6))
    plt.title("Omega vs Moment About x", fontsize=17)
    plt.ylabel("Moment", fontsize=15)
    plt.xlabel("Omega", fontsize=15)
    plt.tick_params(labelsize=13)
    plt.grid()

    mx = [x[0] for x in moments]
    omega = [x[1] for x in variables]
    fz = [f[2] for f in forces]

    scatter = plt.scatter(omega, mx, s=5, c=fz)
    plt.colorbar(scatter, label="$F_z$ (N)")

    plt.legend(fontsize=12)
    plt.show()


def plot_omega_Fz(filename, mass, g):
    # cone = 5
    pitch = -4

    variables, forces, moments = load_data(filename)

    vzs = np.unique([x[3] for x in variables])

    # vzs = [2.5, 2.75, 3, 3.25]

    plt.figure(figsize=(7,5))
    plt.title("Rotation Speed vs Lift Force", fontsize=17)
    plt.ylabel("Force (N)", fontsize=15)
    plt.xlabel("Rotation Speed (rpm)", fontsize=15)
    plt.tick_params(labelsize=13)
    plt.grid()

    plt.text(0.05, 0.9,f"Pitch={pitch}$^{{\circ}}$", transform=plt.gca().transAxes, fontsize=12, bbox={'facecolor': 'white', 'alpha': 0.6})

    # # Initialize lists to store Fz and Vz for each omega
    # fz_data = {vz: [] for vz in vzs}
    # omega_data = {vz: [] for vz in vzs}

    fz = []
    omega = []
    vzs = []

    # Iterate through the data and group by omega
    for i, x in enumerate(variables):

        # if x[0] == cone and x[2] == pitch:
        if x[2] == pitch:
        # if x[0] == cone:
            fz.append(forces[i, 2])  # Fz (lift force)
            omega.append(x[1])         # Vz (descent velocity)
            vzs.append(x[3])

    # fz_arr = np.array(fz_data)
    # omega_arr = np.array(omega_data)

    # for omega in omegas:
    #     vz_squared = np.array(vz_data[omega]) ** 2
    #     fz = np.array(fz_data[omega])
    #     plt.scatter(vz_squared, fz, label=f"omega={omega} rpm")

    scatter = plt.scatter(omega, fz, c=vzs)
    cbar = plt.colorbar(scatter)
    cbar.set_label("$v_z$ (m/s)", size=15)

    plt.savefig(dst+"omegaVsFz.png", bbox_inches='tight', dpi=300)
    plt.close()


def plot_pitch_Fz(filename, mass, g):

    variables, forces, moments = load_data(filename)

    plt.figure(figsize=(10,5))
    plt.title("Pitch vs Force in z", fontsize=17)
    plt.ylabel("Force (N)", fontsize=15)
    plt.xlabel("Pitch Angle (degrees)", fontsize=15)
    plt.tick_params(labelsize=13)
    plt.grid()

    mz = [x[2] for x in moments]
    pitch = [x[2] for x in variables]
    fz = [f[2] for f in forces]

    scatter = plt.scatter(pitch, fz, s=5, c=mz)
    plt.colorbar(scatter, label="$M_z$ (Nm)")

    plt.legend()
    plt.show()


def plot_pitch_Mx(filename, mass, g, I):

    cone = 5

    variables, forces, moments = load_data(filename)

    plt.figure(figsize=(7,5))
    plt.title("Pitch vs Net Moment About x", fontsize=17)
    plt.ylabel("Moment (Nm)", fontsize=15)
    plt.xlabel("Pitch Angle (degrees)", fontsize=15)
    plt.tick_params(labelsize=13)
    plt.grid()

    # mx = [x[0] for x in moments]
    # pitch = [x[2] for x in variables]
    # vz = [x[3] for x in variables]

    vzs = []
    mx = []
    pitch = []
    fz = []
    omegas = []

    vz = 3

    # # Iterate through the data and group by omega
    for i, x in enumerate(variables):

    #     # if x[0] == cone and x[2] == pitch:
    #     if x[0] == cone:
        if x[0] == cone and x[3] == vz:
            # vz.append(x[3])  # Fz (lift force)
            fz.append(forces[i,2])
            omegas.append(x[1])
            pitch.append(x[2])         # Vz (descent velocity)
            m_mass = mass_moment(x, mass, g, I)
            mx.append(moments[i,0]+m_mass[0])

    scatter = plt.scatter(pitch, mx, c=omegas)
    cbar = plt.colorbar(scatter)
    cbar.set_label("Rotation Rate (rpm)", size=15)

    plt.text(0.7, 0.8,f"Cone={cone}$^{{\circ}}$\n$v_z$={vz} m/s", transform=plt.gca().transAxes, fontsize=13, bbox={'facecolor': 'white', 'alpha': 0.6})


    plt.savefig(dst+"pitchVsMx.png", bbox_inches='tight', dpi=300)
    plt.close()


def plot_cone_Mz(filename, mass, g):

    variables, forces, moments = load_data(filename)

    plt.figure(figsize=(7,5))
    plt.title("Cone vs Moment About z", fontsize=17)
    plt.ylabel("Moment (Nm)", fontsize=15)
    plt.xlabel("Cone Angle (degrees)", fontsize=15)
    plt.tick_params(labelsize=13)
    plt.grid()

    mz = [x[2] for x in moments]
    cone = [x[0] for x in variables]
    fz = [f[2] for f in forces]

    scatter = plt.scatter(cone, mz, s=5, c=fz)
    cbar = plt.colorbar(scatter, label="$F_z$ (N)")
    cbar.set_label("$F_z$ (N)", size=15)

    plt.legend()
    plt.show()


def plot_cone_Mx(filename, mass, g):

    variables, forces, moments = load_data(filename)

    plt.figure(figsize=(7,5))
    plt.title("Cone vs Moment About x", fontsize=17)
    plt.ylabel("Moment (Nm)", fontsize=15)
    plt.xlabel("Cone Angle (degrees)", fontsize=15)
    plt.tick_params(labelsize=13)
    plt.grid()

    mx = [x[0] for x in moments]
    cone = [x[0] for x in variables]
    vz = [x[3] for x in variables]

    scatter = plt.scatter(cone, mx, c=vz)
    cbar = plt.colorbar(scatter)
    cbar.set_label("$v_z$ (m/s)", size=15)

    plt.savefig(dst+"coneVsMx.png", bbox_inches='tight', dpi=300)
    plt.close()


def plot_cone_My(filename, mass, g):

    variables, forces, moments = load_data(filename)

    plt.figure(figsize=(7,5))
    plt.title("Cone vs Moment About y", fontsize=17)
    plt.ylabel("Moment (Nm)", fontsize=15)
    plt.xlabel("Cone Angle (degrees)", fontsize=15)
    plt.tick_params(labelsize=13)
    plt.grid()

    mx = [x[1] for x in moments]
    cone = [x[0] for x in variables]
    vz = [x[3] for x in variables]

    scatter = plt.scatter(cone, mx, c=vz)
    cbar = plt.colorbar(scatter)
    cbar.set_label("$v_z$ (m/s)", size=15)

    plt.savefig(dst+"coneVsMy.png", bbox_inches='tight', dpi=300)
    plt.close()


def plot_cone_Mz(filename, mass, g):

    variables, forces, moments = load_data(filename)

    plt.figure(figsize=(7,5))
    plt.title("Cone vs Moment About z", fontsize=17)
    plt.ylabel("Moment (Nm)", fontsize=15)
    plt.xlabel("Cone Angle (degrees)", fontsize=15)
    plt.tick_params(labelsize=13)
    plt.grid()

    mx = [x[2] for x in moments]
    cone = [x[0] for x in variables]
    vz = [x[3] for x in variables]

    scatter = plt.scatter(cone, mx, c=vz)
    cbar = plt.colorbar(scatter)
    cbar.set_label("$v_z$ (m/s)", size=15)

    plt.savefig(dst+"coneVsMz.png", bbox_inches='tight', dpi=300)
    plt.close()


def plot_omega_Mz(filename, mass, g):

    cone = 5

    variables, forces, moments = load_data(filename)

    plt.figure(figsize=(7,5))
    plt.title("Rotation Speed vs Moment About z", fontsize=17)
    plt.ylabel("Moment (Nm)", fontsize=15)
    plt.xlabel("Rotation Speed / $v_z$ (RPM / m/s)", fontsize=15)
    plt.tick_params(labelsize=13)
    plt.grid()

    # mz = [x[2] for x in moments]
    # omega = [x[1] for x in variables]
    # vz = [x[3] for x in variables]


    vz = []
    mz = []
    omega = []
    # pitch = []

    # vz = 3
    pitch = -4

    # # Iterate through the data and group by omega
    for i, x in enumerate(variables):

        # if x[2] == pitch:
        # if x[0] == cone:
        if x[2] == pitch and x[0] == cone:  # and x[3] < 3.5:
        
            # pitch.append(x[2])
            vz.append(x[3])  # Fz (lift force)
            omega.append(x[1])         # Vz (descent velocity)
            mz.append(moments[i,2])

    omega = np.array(omega)
    vz = np.array(vz)

    scatter = plt.scatter(omega / vz, mz, c=vz)
    cbar = plt.colorbar(scatter)
    cbar.set_label("$v_z$ (m/s)", size=15)

    print("max", omega[np.argmax(mz)])

    plt.text(0.7, 0.8,f"Cone={cone}\nPitch={pitch}", transform=plt.gca().transAxes, fontsize=13, bbox={'facecolor': 'white', 'alpha': 0.6})


    plt.savefig(dst+"omegaVsMz.png", bbox_inches='tight', dpi=300)
    plt.close()


def plot_omega_Mx(filename, mass, g):

    cone = 5

    variables, forces, moments = load_data(filename)

    plt.figure(figsize=(7,5))
    plt.title("Rotation Speed vs Moment About x", fontsize=17)
    plt.ylabel("Moment (Nm)", fontsize=15)
    plt.xlabel("Rotation Speed / $v_z$ (RPM / m/s)", fontsize=15)
    plt.tick_params(labelsize=13)
    plt.grid()

    # mz = [x[2] for x in moments]
    # omega = [x[1] for x in variables]
    # vz = [x[3] for x in variables]


    vz = []
    mz = []
    omega = []
    # pitch = []
    cone = []

    # vz = 3
    pitch = -4

    # # Iterate through the data and group by omega
    for i, x in enumerate(variables):

        # if x[2] == pitch:
        #if x[0] == cone:
        if x[2] == pitch and x[3] == 3.25:
        
            # pitch.append(x[2])
            vz.append(x[3])  # Fz (lift force)
            omega.append(x[1])         # Vz (descent velocity)
            mz.append(moments[i,0])
            cone.append(x[0])

    omega = np.array(omega)
    vz = np.array(vz)

    scatter = plt.scatter(omega, mz, c=vz)
    cbar = plt.colorbar(scatter)
    cbar.set_label("$v_z$ (m/s)", size=15)

    print("max", omega[np.argmax(mz)])

    plt.text(0.7, 0.8,f"Cone={cone}\nPitch={pitch}", transform=plt.gca().transAxes, fontsize=13, bbox={'facecolor': 'white', 'alpha': 0.6})


    plt.savefig(dst+"omegaVsMx.png", bbox_inches='tight', dpi=300)
    plt.close()



def plot_cone_Fz(filename, mass, g):

    cone = 5

    variables, forces, moments = load_data(filename)

    plt.figure(figsize=(7,5))
    plt.title("Cone Angle vs Lift Force", fontsize=17)
    plt.ylabel("Force (N)", fontsize=15)
    plt.xlabel("Pitch Angle (degrees)", fontsize=15)
    plt.tick_params(labelsize=13)
    plt.grid()

    # mz = [x[2] for x in moments]
    # omega = [x[1] for x in variables]
    # vz = [x[3] for x in variables]


    vz = []
    mz = []
    omega = []
    pitch = []
    cone = []
    fz = []

    vz = 3

    # # Iterate through the data and group by omega
    for i, x in enumerate(variables):

        # if x[0] == cone and x[2] == pitch:
        # if x[0] == cone:
        if x[3] == vz:
            cone.append(x[0])
            fz.append(forces[i,2])  # Fz (lift force)
            # omega.append(x[1])         # Vz (descent velocity)
            mz.append(moments[i,2])
            pitch.append(x[2])

    scatter = plt.scatter(pitch, fz, c=cone)
    cbar = plt.colorbar(scatter)
    cbar.set_label("cone (degrees)", size=15)

    plt.savefig(dst+"coneVsFz.png", bbox_inches='tight', dpi=300)
    plt.close()


def plot_Mnet_omega(filename, mass, g, I):
    
    cone = 5

    variables, forces, moments = load_data(filename)

    plt.figure(figsize=(7,5))
    plt.title("Cone Angle vs Lift Force", fontsize=17)
    plt.ylabel("Force (N)", fontsize=15)
    plt.xlabel("Pitch Angle (degrees)", fontsize=15)
    plt.tick_params(labelsize=13)
    plt.grid()

    # mz = [x[2] for x in moments]
    # omega = [x[1] for x in variables]
    # vz = [x[3] for x in variables]


    vz = []
    mz = []
    omega = []
    pitches = []
    cone = []
    fz = []
    netM = []


    cone = 5
    pitch = -4

    # # Iterate through the data and group by omega
    for i, x in enumerate(variables):

        if x[0] == cone and x[2] == pitch:

            M_inert = mass_moment(x)

            pitches.append(x[2])

            fz.append(np.abs(forces[i,2] + M_inert[3]))  # Fz (lift force)
            omega.append(x[1])         # Vz (descent velocity)
            M = np.linalg.norm( moments[i,:] + M_inert[:3] )

            netM.append(M)
            #pitch.append(x[2])

    scatter = plt.scatter(omega, netM, c=fz)
    cbar = plt.colorbar(scatter)
    cbar.set_label("|Net $F_z$| (N)", size=15)

    plt.set_cmap('coolwarm')

    plt.savefig(dst+"omegaVsNetM.png", bbox_inches='tight', dpi=300)
    plt.close()



def plot_omega_My(filename, mass, g, I):

    cone = 5
    pitch = -4

    variables, forces, moments = load_data(filename)

    plt.figure(figsize=(7,5))
    plt.title("Rotation Speed vs Net Moment About y", fontsize=17)
    plt.ylabel("Moment (Nm)", fontsize=15)
    plt.xlabel("Rotation Speed (rpm)", fontsize=15)
    plt.tick_params(labelsize=13)
    plt.grid()

    # mx = [x[0] for x in moments]
    # pitch = [x[2] for x in variables]
    # vz = [x[3] for x in variables]

    vzs = []
    my = []

    fz = []
    omegas = []

    # vz = 3

    # # Iterate through the data and group by omega
    for i, x in enumerate(variables):

    #     # if x[0] == cone and x[2] == pitch:
    #     if x[0] == cone:
        if x[0] == cone and x[2] == pitch:
            vzs.append(x[3])  # Fz (lift force)
            fz.append(forces[i,2])
            omegas.append(x[1]*x[1])

            m_mass = mass_moment(x, mass, g, I)
            my.append(moments[i,1]+m_mass[1])

    scatter = plt.scatter(omegas, my, c=vzs)
    cbar = plt.colorbar(scatter)
    cbar.set_label("$v_z$ (m/s)", size=15)

    plt.text(0.05, 0.8,f"Cone={cone}$^{{\circ}}$\nPitch={pitch}$^{{\circ}}$", transform=plt.gca().transAxes, fontsize=13, bbox={'facecolor': 'white', 'alpha': 0.6})

    plt.savefig(dst+"omegaVsMy.png", bbox_inches='tight', dpi=300)
    plt.close()


def plot_omega_My2(filename, cone, pitch, vzs, I):
    """
    """

    # Load data from the file
    variables, forces, moments = load_data(filename)

    # Initialize lists to store Fz and Vz for each omega

    omega_data = {vz: [] for vz in vzs}
    my_data = {vz: [] for vz in vzs}

    # Iterate through the data and group by omega
    for i, x in enumerate(variables):
        for vz in vzs:
            if x[0] == cone and x[2] == pitch and x[3] == vz: # and x[1] >= 200:
                m_mass = mass_moment(x, mass, g, I)
                my_data[vz].append(moments[i, 1]+m_mass[1])  # Fz (lift force)
                omega_data[vz].append(x[1])         # Vz (descent velocity)

    # Plot Fz vs Vz^2 for each omega
    plt.figure(figsize=(7,5))
    plt.title("Moment About y Scaling with Rotation Speed", fontsize=17)
    plt.ylabel("Moment (Nm)", fontsize=15)
    plt.xlabel("Rotation Speed Squared (rpm$^2$)", fontsize=15)
    plt.tick_params(labelsize=13)

    for vz in vzs:
        omega = np.array(omega_data[vz])
        omega_squared = np.array(omega) ** 2
        omega_squared = omega_squared/vz**2  # *np.sin(cone*np.pi/180) - np.cos(cone*np.pi/180)*vz**2
        my = np.array(my_data[vz])

        # Fit a linear regression line
        if len(omega_squared) > 1:  # Ensure there are enough points to fit
            slope, intercept, r_value, _, _ = linregress(omega_squared, my)
            my_fit = slope * omega_squared + intercept
            plt.plot(omega_squared, my_fit, ls="-", label=f"{vz} m/s\n$R^2$={r_value**2:.3f}")
        plt.scatter(omega_squared, my)


    plt.text(0.7, 0.05,f"Cone={cone}$^{{\circ}}$\nPitch={pitch}$^{{\circ}}$", transform=plt.gca().transAxes, fontsize=12, bbox={'facecolor': 'white', 'alpha': 0.6})

    plt.grid()
    plt.legend(fontsize=13)
    plt.savefig(dst+"omega-my.png", bbox_inches="tight", dpi=300)
    plt.close()




if __name__ == '__main__':

    filename = sys.argv[1]

    dst = "test_figures/"
    
    cone = 5
    pitch = 0
    omegas = [240,450,600]
    vzs=[2,2.5,3]

    inertia=[0.000013713365031784028, 0.0000898853777707176, 0.00010289524930712706]

    # loading inertia about principle axes into 3x3 array
    I = np.zeros(shape=(3,3))

    I[0,0] = inertia[0]
    I[1,1] = inertia[1]
    I[2,2] = inertia[2]

    # plot_fz_vz(filename, cone, pitch, omegas)

    mass = 0.0281
    g = 9.8067

    # plot_fz_vz(filename, cone, pitch, omegas)
    # plot_pitch_Mz(filename, mass, g)
    # plot_omega_Mx(filename, mass, g)
    # # # plot_pitch_omega(filename, mass, g)
    # plot_omega_Fz(filename, mass, g)
    # # plot_pitch_Fz(filename, mass, g)
    # # # plot_cone_Mz(filename, mass, g)
    # plot_pitch_Mx(filename, mass, g, I)
    # # plot_cone_Mx(filename, mass, g)
    # # plot_cone_My(filename, mass, g)
    # # plot_cone_Mz(filename, mass, g)
    # plot_omega_Mz(filename, mass, g)
    # # plot_omega_Mx(filename, mass, g)

    # plot_cone_Fz(filename, mass, g)
    # plot_Mnet_omega(filename, mass, g, I)
    # drop_test_result()
    # plot_omega_My(filename, mass, g, I)
    plot_omega_My2(filename, cone, pitch, vzs, I)