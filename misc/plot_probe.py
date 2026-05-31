"""
plotting theta - omega
"""

import numpy as np
import matplotlib.pyplot as plt
import argparse
import sys


def load_data(filepath):
    """
    input: path to data file\n
    output: x array, y arrays for number of columns in data
    """

    data = np.loadtxt(filepath)
    
    # loads first column x data, and all other columns as y data(s?)
    x = data[:,0]
    ys = data[:,1:]

    return x, ys


def plot_small_angle(filepath):
    """
    input: path to data file\n
    output: save plot of mass vs eigenfrequencies
    """

    # theta and the numerical and analytical solutions
    x, ys = load_data(filepath)

    root = x[0]*180/np.pi

    x = x[1:]
    ys = ys[1:,:]

    x_deg = np.rad2deg(x)


    plt.figure(figsize=(10,6))

    plt.title(r"f($\theta$) - $\theta$", fontsize=18)
    plt.ylabel(r"f($\theta$)", fontsize=16)
    plt.xlabel(r"$\theta$ (deg)", fontsize=16)

    plt.plot(x_deg, ys[:,0], label=f"Numerical")
    plt.plot(x_deg, ys[:,2], label=f"Analytical - 2nd order")
    plt.plot(x_deg, ys[:,1], label=f"6th order")
    

    plt.axvline(root, ls='--', color='r', label=fr"$\theta_0$ = {root:.3f}")


    plt.grid()
    plt.legend(loc='upper left',fontsize='large')

    # increase tick size
    plt.xticks(fontsize=14)
    plt.yticks(fontsize=14)

    plt.savefig("results/theta-paper.png")
    plt.close()

    # exit immediately (ensure no long-running code after save)
    sys.stdout.flush()
    sys.exit(0)


if __name__ == '__main__':
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
        plot_small_angle(fpath)

