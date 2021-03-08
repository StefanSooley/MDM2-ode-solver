#!/usr/bin/env python3
import numpy as np
from scipy.integrate import solve_ivp
from plot import *
from methods import *


def main():
    args = [m1, m2, l, mu, g, r]
    initial_conds = [x0, phi0, v0, w0]
    print("Integrating...")

    sol = solve_ivp(diff_func_L5, (t_0, t_end), initial_conds, method='RK45', args=([args]), dense_output=True, max_step = 0.0001)

    print("Integration success = {}, {}".format(sol.success, sol.message))

    plot_solution_panels(sol)
    plot_spiral(sol, l, r)
    plot_colour_spiral(sol, l, r)


if __name__ == '__main__':
    # Initial Conditions
    x0 = 0
    phi0 = 0
    v0 = 0
    w0 = 0
    # Constant values
    m1 = 250 # Heavier mass
    m2 = 50  # Lighter mass
    l = 20
    mu = 0.5
    g = 9.81
    r = 1  # Radius of the "pulley"
    # Integration region
    t_0 = 0
    t_end = 15

    main()
