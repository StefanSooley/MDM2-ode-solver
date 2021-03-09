#!/usr/bin/env python3
import matplotlib.pyplot as plt
from matplotlib import colors, cm
import numpy as np
from file_manage import *

plt.rcParams["mathtext.fontset"] = "cm"
plt.rcParams["font.family"] = "Times New Roman"


def plot_solution_panels(sol):
    print("Plotting solution panels...")
    fig, axs = plt.subplots(2, 2)

    axs[0, 0].plot(sol.t, sol.y[0])
    axs[0, 0].set_title('$x$ against $t$')
    axs[0, 0].set_ylabel('$x$ $(m)$')
    axs[0, 0].grid()

    axs[0, 1].plot(sol.t, sol.y[1], 'tab:orange')
    axs[0, 1].set_title('$\phi$ against $t$')
    axs[0, 1].set_ylabel('$\phi$ $(rads)$')
    axs[0, 1].grid()

    axs[1, 0].plot(sol.t, sol.y[2], 'tab:green')
    axs[1, 0].set_title('$v$ against $t$')
    axs[1, 0].set_ylabel('$v$ $(ms^{-1})$')
    axs[1, 0].grid()

    axs[1, 1].plot(sol.t, sol.y[3], 'tab:red')
    axs[1, 1].set_title('$\omega$ against $t$')
    axs[1, 1].set_ylabel('$\omega$ $(radss^{-1})$')
    axs[1, 1].grid()

    fig.tight_layout()

    plt.show()


def plot_spiral(sol, l, r):
    print("Plotting m2 spiral...")
    ax = plt.subplot(111, polar=True)

    phi_arr = sol.y[1]
    x_arr = sol.y[0]
    LHS_arr = np.zeros(len(phi_arr))

    for idx in range(len(phi_arr)):
        LHS = l - phi_arr[idx] * r - x_arr[idx]
        if LHS > 0:
            LHS_arr[idx] = LHS
        else:
            LHS_arr[idx] = 0

    ax.plot(phi_arr, LHS_arr, 'r')
    ax.set_theta_offset(np.pi / 2)
    ax.set_xticklabels([r'$\frac{3\pi}{2}$', r'$\frac{7\pi}{4}$', '0', r'$\frac{\pi}{4}$', r'$\frac{\pi}{2}$',
                        r'$\frac{3\pi}{4}$', r'$\pi$', r'$\frac{5\pi}{4}$'])

    ax.set_title("Trajectory of $m$")

    plt.show()


def plot_colour_spiral(sol, l, r):
    print("Plotting m2 coloured spiral...")
    x_arr = sol.y[0]
    phi_arr = sol.y[1]
    w_arr = sol.y[3]
    log_w_arr = np.zeros(len(phi_arr))
    lhs_arr = np.zeros(len(phi_arr))
    colour_range = (-5, 2)
    # Need to calculate the length of the LHS
    for idx in range(len(phi_arr)):
        LHS = l - phi_arr[idx] * r - x_arr[idx]
        if LHS > 0:
            lhs_arr[idx] = LHS
        else:
            lhs_arr[idx] = 0

    # To make the colours on the line less extreme, they are logged as the range of values is very high
    max_val = max(w_arr)
    for idx, i in enumerate(w_arr):
        if i > 0:
            log_w_arr[idx] = (np.log(i))
        else:
            log_w_arr[idx] = (np.log(max_val))

    colormap = plt.get_cmap('inferno')
    norm = colors.Normalize(*colour_range)
    ax = plt.subplot(1, 1, 1, polar=True)

    # Colour bar
    cbar = plt.colorbar(cm.ScalarMappable(norm=norm, cmap=colormap), ax=ax)
    ticklabels = ["min $\omega$", "max $\omega$"]
    cbar.set_ticks(np.linspace(*colour_range, len(ticklabels)))
    cbar.set_ticklabels(ticklabels)

    # The actual line, its a scatter graph where each point is the colour based off the w at that time.
    ax.scatter(phi_arr, lhs_arr, c=log_w_arr, s=10, cmap=colormap, norm=norm, linewidths=0)
    ax.set_theta_offset(np.pi / 2)
    ax.set_xticklabels([r'$\frac{3\pi}{2}$', r'$\frac{7\pi}{4}$', '0', r'$\frac{\pi}{4}$', r'$\frac{\pi}{2}$',
                        r'$\frac{3\pi}{4}$', r'$\pi$', r'$\frac{5\pi}{4}$'])
    ax.set_title("Trajectory of $m$")
    ax.grid(True)

    plt.show()


def plot_multi_sol_panels(any_sol_type, legend_labels):
    """
            Plots each variable in the solution over t, from different solutions of the ODE on the same graph.

            Arguments:
                any_sol_type: Either the filename of the array of solutions or the solution array itself.

                legend_labels: An array of the labels of the lines to be shown in the legend. There must be the same
                number of labels as there are solutions:
                            legend_labels = ["name1", "name2"]

            """

    if type(any_sol_type) is str:
        sol_arr = open_file(any_sol_type)
    else:
        sol_arr = any_sol_type
    if len(legend_labels) != len(sol_arr):
        raise Exception("There arent the same number of labels as solutions")

    print("Plotting solution panels...")
    fig, axs = plt.subplots(2, 2)

    axs[0, 0].set_title('$x$ against $t$')
    axs[0, 0].set_ylabel('$x$ $(m)$')
    axs[0, 0].grid()
    axs[0, 1].set_title('$\phi$ against $t$')
    axs[0, 1].set_ylabel('$\phi$ $(rads)$')
    axs[0, 1].grid()
    axs[1, 0].set_title('$v$ against $t$')
    axs[1, 0].set_ylabel('$v$ $(ms^{-1})$')
    axs[1, 0].grid()
    axs[1, 1].set_title('$\omega$ against $t$')
    axs[1, 1].set_ylabel('$\omega$ $(radss^{-1})$')
    axs[1, 1].grid()

    for idx, sol in enumerate(sol_arr):
        axs[0, 0].plot(sol.t, sol.y[0], label=legend_labels[idx])
        axs[0, 1].plot(sol.t, sol.y[1], label=legend_labels[idx])
        axs[1, 0].plot(sol.t, sol.y[2], label=legend_labels[idx])
        axs[1, 1].plot(sol.t, sol.y[3], label=legend_labels[idx])
    axs[0, 0].legend()
    axs[0, 1].legend()
    axs[1, 0].legend()
    axs[1, 1].legend()
    fig.tight_layout()
    plt.show()


def plot_multi_spiral(any_sol_type, legend_labels, args_arr):
    """
            Plots multiple trajectories of m2, from different solutions of the ODE.

            Arguments:
                any_sol_type: Either the filename of the array of solutions or the solution array itself.

                legend_labels: An array of the labels of the lines to be shown in the legend. There must be the same
                number of labels as there are solutions:
                            legend_labels = ["name1", "name2"]

                args_arr: The array of the arguments used to find the solutions plotted. Required for calculating the
                LHS of the length.

            """


    if type(any_sol_type) is str:
        sol_arr = open_file(any_sol_type)
    else:
        sol_arr = any_sol_type
    if len(legend_labels) != len(sol_arr):
        raise Exception("There arent the same number of labels as solutions")

    print("Plotting m2 spiral...")
    ax = plt.subplot(111, polar=True)
    for idx, sol in enumerate(sol_arr):
        name = legend_labels[idx]
        phi_arr = sol.y[1]
        x_arr = sol.y[0]
        LHS_arr = np.zeros(len(phi_arr))
        l = args_arr[idx][2]
        r = args_arr[idx][5]

        for idx in range(len(phi_arr)):
            LHS = l - phi_arr[idx] * r - x_arr[idx]
            if LHS>0:
                LHS_arr[idx] = LHS
            else:
                LHS_arr[idx] = 0

        ax.plot(phi_arr, LHS_arr, label=name)
    ax.set_theta_offset(np.pi / 2)
    ax.set_xticklabels([r'$\frac{3\pi}{2}$', r'$\frac{7\pi}{4}$', '0', r'$\frac{\pi}{4}$', r'$\frac{\pi}{2}$',
                        r'$\frac{3\pi}{4}$', r'$\pi$', r'$\frac{5\pi}{4}$'])

    ax.legend()
    ax.set_title("Trajectory of $m$")

    plt.show()
