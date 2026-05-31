"""
for plotting the trajectories of venus probes from the balloon
"""


import numpy as np
import matplotlib.pyplot as plt
from matplotlib.animation import FuncAnimation
from tqdm import tqdm
import sys


def load_data(filename):
    """
    """

    global n_probes

    with open(filename) as f:
        topline = f.readline().strip("\n")
        
        topline = topline.split()

        n_probes = int(topline[-1])

    data = np.loadtxt(filename, comments='#')

    time = data[:,0]
    data = data[:,1:]

    print(time[0:10])

    return time, data


def animate(filename):
    """"""

    time, data = load_data(filename)

    timestep = time[1] - time[0]
    timesteps = len(time)

    x = data[:,0:timesteps:3]
    y = data[:,1:timesteps:3]
    z = data[:,2:timesteps:3]

    fig = plt.figure(figsize=(12,8))
    ax = fig.add_subplot(projection='3d')
    fig.suptitle("Probe Trajectories")

    # text object to tell the time and energies on the video
    fig_text = fig.text(0.05, 0.9, "", fontsize=11)

    # set axis limits once
    ax.set_xlim(x.min(), x.max())
    ax.set_ylim(y.min(), y.max())
    ax.set_zlim(z.min(), z.max())

    # PROGRESS BAR FOR FUN
    pbar = tqdm(total=timesteps, desc="Rendering frames", ncols=100)

    # create a line (trail) + point for each mass
    trails = []
    points = []
    for i in range(n_probes):
        line, = ax.plot([], [], [], "-", linewidth=1)   # trail
        pt,   = ax.plot([], [], [], "o")                # particle
        trails.append(line)
        points.append(pt)

    def update(frame):

        for i in range(n_probes):
            x_trail = x[:frame+1, i]
            y_trail = y[:frame+1, i]
            z_trail = z[:frame+1, i]

            trails[i].set_data(x_trail, y_trail)
            trails[i].set_3d_properties(z_trail)        # full trail

            points[i].set_data([x_trail[-1]], [y_trail[-1]])   # last point only
            points[i].set_3d_properties([z_trail[-1]])

        fig_text.set_text(
            f"t = {time[frame]:.2f} s"
        )

        pbar.update(1)
        return trails + points
    
    filename = filename.split('.')
    filename = str(filename[0])
    
    ani = FuncAnimation(fig, update, frames=timesteps, interval=40, blit=True)
    ani.save(f"video_{filename}.mp4", writer="ffmpeg", fps=60)
    plt.close()


def plot_trajectory(filename):
    
    time, data = load_data(filename)

    timestep = time[1] - time[0]
    timesteps = len(time)

    x = data[:,0:timesteps:3]
    y = data[:,1:timesteps:3]
    z = data[:,2:timesteps:3]

    x = np.array(x) / 1000
    y = np.array(y) / 1000

    z = np.array(z)
    z = z / 1000

    fig = plt.figure(figsize=(12,8))
    ax = fig.add_subplot(projection='3d')
    fig.suptitle("Probe Trajectories")

    fig.tight_layout()

    # set axis limits once
    ax.set_xlim(x.min(), x.max())
    ax.set_ylim(y.min(), y.max())
    ax.set_zlim(z.min(), z.max())

    for i in range(n_probes):
        ax.scatter(x[:,i], y[:,i], z[:,i])
        ax.plot(x[:,i], y[:,i], z[:,i])

    plt.savefig("probePath3D.png", dpi=300)
    plt.close()

    fig, axs = plt.subplots(1,2)

    for i in range(n_probes):
        axs[0].plot(x[:,i],z[:,i],label=f'xz_{i}')
        axs[1].plot(y[:,i],z[:,i],label=f'yz_{i}')

    # axs[0].legend()
    # axs[1].legend()

    plt.savefig("probePath2D.png", dpi=300)
    plt.close()


    plt.title("Probe Distance From Balloon", fontsize=16)
    
    rho = np.zeros((timesteps,n_probes))

    for i in range(timesteps):
        rho[i,:] = np.sqrt(x[i,:]**2 + y[i,:]**2)

    for i in range(n_probes):
        plt.plot(rho[:,i], z[:,i], label=f'Probe-{i}')

    plt.xlabel("Radial distance from deployment site (km)", fontsize=14)
    plt.ylabel("Altitude (km)", fontsize=14)
    plt.tick_params(labelsize=12)

    # plt.legend()
    plt.grid()
    plt.savefig("probeDistance.png", dpi=300)
    plt.close()


if __name__ == '__main__':

    file = sys.argv[1]

    plot_trajectory(file)

    # animate(file)
