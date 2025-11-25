"""
rotation of an object

for the plotting of the results of the rigid body integration in rk4.f90
"""


import trimesh
from mpl_toolkits.mplot3d.art3d import Poly3DCollection
import numpy as np
import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d import Axes3D
from matplotlib.animation import FuncAnimation
import argparse
from tqdm import tqdm


def load_data(filename):
    """input: filename\n
    output: omega in x, y, z."""

    omegas = np.loadtxt(filename)

    # wx = data[:,0]
    # wy = data[:,1]
    # wz = data[:,2]

    # has the capability to receive more columns (torque, etc.)

    # return wx, wy, wz
    return omegas

# def rotation_matrix(omega, t):
#     """return combined rotation matrix for angular speeds (omega_x, omega_y, omega_z) and time t."""
#     wx, wy, wz = omega

#     Rx = np.array([[1,            0,             0],
#                    [0, np.cos(wx*t), -np.sin(wx*t)],
#                    [0, np.sin(wx*t),  np.cos(wx*t)]])

#     Ry = np.array([[ np.cos(wy*t), 0, np.sin(wy*t)],
#                    [0            ,1            , 0],
#                    [-np.sin(wy*t), 0, np.cos(wy*t)]])

#     Rz = np.array([[np.cos(wz*t), -np.sin(wz*t), 0],
#                    [np.sin(wz*t),  np.cos(wz*t), 0],
#                    [0,           0,              1]])

#     # R is the rotation matrix, matmul product of Rx,y,z (in reverse order for euler)
#     R = Rz @ Ry @ Rx

#     return R

def rotation_matrix(omega, dt):
    """Return rotation matrix for angular velocity vector omega over dt seconds."""
    theta = np.linalg.norm(omega) * dt
    if theta == 0:
        return np.eye(3)
    k = omega / np.linalg.norm(omega)
    K = np.array([[0, -k[2], k[1]],
                  [k[2], 0, -k[0]],
                  [-k[1], k[0], 0]])
    R = np.eye(3) + np.sin(theta) * K + (1 - np.cos(theta)) * (K @ K)
    return R


def init_rod(L=1.0):

    rod = np.array([[0, -L/2, 0], [0, +L/2, 0]])

    R_total = np.eye(3)  # to account for cumulative changes

    return rod, R_total


def animate(filename, stl_path='probe_designs/base_design.stl'):

    omega_data = load_data(filename)

    # number of frames = number of iterations in rk4.f90
    n_frames = len(omega_data[:,0])

    t = omega_data[:,0]


    # this bit will probably go in favour of reading in an stl
    # L = 1.5
    # rod, R_total = init_rod(L)

    # load stl
    mesh = trimesh.load(stl_path)
    vertices = mesh.vertices
    faces = mesh.faces

    # mesh.apply_translation(-mesh.centroid)

    # init rotation matrix
    R_total = np.eye(3)

    fig = plt.figure()
    ax = fig.add_subplot(111, projection='3d')

    poly = Poly3DCollection(vertices[faces], alpha=0.8)
    ax.add_collection3d(poly)

    # scales
    scale = vertices.flatten()
    # ax.auto_scale_xyz(scale, scale, scale)

    # to handle probe rotating out for bounds
    L = np.max(np.abs(vertices)) * 1.6

    ax.set_xlim(-L, L)
    ax.set_ylim(-L, L)
    ax.set_zlim(-L, L)
    ax.set_box_aspect([1, 1, 1])
    ax.set_xlabel('X')
    ax.set_ylabel('Y')
    ax.set_zlabel('Z')

    # AXIS TO ROTATE WITH THE PROBE TO SEE ORIENTATION
    triad_len = L * 0.8  # scale relative to model size
    # ax.plot([0, triad_len], [0, 0], [0, 0], color='r', lw=2, label='X-axis')
    # ax.plot([0, 0], [0, triad_len], [0, 0], color='g', lw=2, label='Y-axis')
    # ax.plot([0, 0], [0, 0], [0, triad_len], color='b', lw=2, label='Z-axis')
    # ax.legend()
    
    body_axes = [
    ax.plot([0, triad_len], [0, 0], [0, 0], color='r', lw=2)[0],  # X
    ax.plot([0, 0], [0, triad_len], [0, 0], color='g', lw=2)[0],  # Y
    ax.plot([0, 0], [0, 0], [0, triad_len], color='b', lw=2)[0],  # Z
    ]
    # (line,) = ax.plot([], [], [], 'r-', lw=3)

    # text object to tell the time on the video
    time_text = fig.text(0.05, 0.8, "", fontsize=12)

    # PROGRESS BAR FOR FUN
    pbar = tqdm(total=n_frames, desc="Rendering frames", ncols=80)

    def update(i):
        nonlocal R_total
        omega = omega_data[i,1:]
        if i == 0:
            pbar.update(1)
            return (poly,)
        dt = t[i] - t[i - 1]
        R_step = rotation_matrix(omega, dt)
        R_total = R_step @ R_total
        rotated = vertices @ R_total.T
        # line.set_data(rod_rot[:,0], rod_rot[:,1])
        # line.set_3d_properties(rod_rot[:,2])
        # return (line,)
        poly.set_verts(rotated[faces])

        origin = np.array([[0, 0, 0]])
        axes_dirs = np.eye(3) * triad_len
        rotated_axes = axes_dirs @ R_total.T

        for j in range(3):
            xdata = [origin[0, 0], rotated_axes[j, 0]]
            ydata = [origin[0, 1], rotated_axes[j, 1]]
            zdata = [origin[0, 2], rotated_axes[j, 2]]
            body_axes[j].set_data(xdata, ydata)
            body_axes[j].set_3d_properties(zdata)

        time_text.set_text(
            f"t = {t[i]:.3f} s\n"
            fr"$\omega$ = [{{{omega[0]:.2f}}}, {{{omega[1]:.2f}}}, {{{omega[2]:.2f}}}] rad/s"
        )

        pbar.update(1)
        return (poly,)

    # do the animation
    ani = FuncAnimation(fig, update, frames=n_frames, interval=40, blit=True)
    ani.save("results/rotation.mp4", fps=30, dpi=150, codec="libx264")

    pbar.close()


def animate_multi(filename, stl_file='base_design.stl'):
    omega_data = load_data(filename)

    # number of frames = number of iterations in rk4.f90
    n_frames = len(omega_data[:,0])

    t = omega_data[:,0]
    omega_data = omega_data[:,1:]

    mesh = trimesh.load(stl_file)
    vertices = np.array(mesh.vertices)
    faces = np.array(mesh.faces)

    L = 1.2 * np.max(np.abs(vertices))

    # ----------- Create 3 subplots -----------
    fig = plt.figure(figsize=(12, 4))

    # text object to tell the time on the video
    time_text = fig.text(0.05, 0.8, "", fontsize=12)


    ax1 = fig.add_subplot(131, projection='3d')
    ax2 = fig.add_subplot(132, projection='3d')
    ax3 = fig.add_subplot(133, projection='3d')

    for ax in (ax1, ax2, ax3):
        ax.set_xlim(-L, L)
        ax.set_ylim(-L, L)
        ax.set_zlim(-L, L)

    ax1.set_title("All over view")
    ax2.set_title("Top (Z)")
    ax3.set_title("Side (Y)")

    ax1.view_init(elev=30, azim=45)   # 3D angle  
    ax2.view_init(elev=90, azim=0)    # top-down  
    ax3.view_init(elev=0, azim=0)     # side  

    # Create trisurf placeholders for each plot
    plot1 = [ax1.plot_trisurf(vertices[:,0], vertices[:,1], vertices[:,2], triangles=faces, color='orange', alpha=0.9)]
    plot2 = [ax2.plot_trisurf(vertices[:,0], vertices[:,1], vertices[:,2], triangles=faces, color='orange', alpha=0.9)]
    plot3 = [ax3.plot_trisurf(vertices[:,0], vertices[:,1], vertices[:,2], triangles=faces, color='orange', alpha=0.9)]

    R_total = np.eye(3)

    # PROGRESS BAR FOR FUN
    pbar = tqdm(total=n_frames, desc="Rendering frames", ncols=80)

    def update(i):
        nonlocal R_total
        omega = omega_data[i]
        if i > 0:
            dt = t[i] - t[i-1]
            omega = omega_data[i]
            R_total = rotation_matrix(omega, dt) @ R_total

        verts_rot = vertices @ R_total.T

        # Remove old meshes
        for pl in (plot1, plot2, plot3):
            pl[0].remove()

        # Re-plot updated mesh
        plot1[0] = ax1.plot_trisurf(verts_rot[:,0], verts_rot[:,1], verts_rot[:,2], triangles=faces, color='blue', alpha=0.9)
        plot2[0] = ax2.plot_trisurf(verts_rot[:,0], verts_rot[:,1], verts_rot[:,2], triangles=faces, color='blue', alpha=0.9)
        plot3[0] = ax3.plot_trisurf(verts_rot[:,0], verts_rot[:,1], verts_rot[:,2], triangles=faces, color='blue', alpha=0.9)

        pbar.update(1)

        time_text.set_text(
            f"t = {t[i]:.3f} s\n"
            fr"$\omega$ = [{{{omega[0]:.2f}}}, {{{omega[1]:.2f}}}, {{{omega[2]:.2f}}}] rad/s"
        )

        return plot1 + plot2 + plot3

    ani = FuncAnimation(fig, update, frames=n_frames, interval=40, blit=False)
    ani.save("results/multi_view_rotation.mp4", fps=30, dpi=150)

    pbar.close()


if __name__ == "__main__":

    parser = argparse.ArgumentParser(description='''
                                     plotting the data output from coupled
                                     harmonic oscillator code''')
    parser.add_argument('filepath', type=str, help='path to data file')

    parser.add_argument('option', type=int, help='option to plot',
                        nargs='?', default="0")
    args = parser.parse_args()

    fpath = args.filepath

    option = args.option

    if option==0:
        animate(fpath)
    if option==1:
        animate_multi(fpath)
