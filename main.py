#!/usr/bin/env python3
import numpy as np
from scipy.integrate import *
from plot import *
from methods import *

# Initial Conditions
x0 = 0.001
v0 = 0
phi0 = 0
w0 = 0
# Constant values
m1 = 10
m2 = 1
l = 20
mu = 0.6
g = 9.81
# Integration region
t_0 = 0
t_end = 10


def main():

    args = [m1, m2, l, mu, g]
    initial_conds = [x0, phi0, v0, w0]

    sol = solve_ivp(diff_func_L2, (t_0, t_end), initial_conds, method='BDF', args=([args]), dense_output=True)

    print("Success = {}, {}".format(sol.success,sol.message))

    plot_solution_panels(sol)


if __name__ == '__main__':
    main()
