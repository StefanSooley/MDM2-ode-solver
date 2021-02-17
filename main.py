#!/usr/bin/env python3
import numpy as np
from scipy.integrate import *
from plot import *


def diff_func(t, vec, args):
    """
        Defines the differential equations for the pulley system.

        Arguments:
            vec :  vector of the state variables:
                      w = [x,phi,v,w]

            t :  time

            args : vector of the constant values:
                      p = [m1, m2, l, mu, g]

    """

    x, phi, v, w = vec
    m1, m2, l, mu, g = args

    rhs = (m1 * g) / np.exp(mu * (phi + (np.pi / 2)))

    delta_x = l - x
    print("delta_x = ", delta_x)
    a = -m1 - m2
    b = -delta_x * m2 * np.sin(phi)
    c = m1 * g - m2 * v * w * np.sin(phi) - m2 * (w ** 2) * delta_x - m2 * g * np.sin(phi) + m2 * w * v * np.sin(
        phi) - m2 * (w ** 2) * delta_x * np.cos(phi) - rhs

    d = -m2 * delta_x * np.sin(phi)
    e = -m2 * (delta_x ** 2)
    f = m2 * g * delta_x * np.cos(phi) + m2 * (v ** 2) * np.sin(phi) + 2 * m2 * v * w * delta_x - rhs

    d_dt = [
        v,
        w,
        (f * b - c * e) / (a * e - d * b),
        (f * a - c * d) / (b * d - e * a)
    ]

    return d_dt


def main():
    # Initial Conditions
    x0 = 0
    v0 = 0
    phi0 = (np.pi) / 2
    w0 = 0
    # Constant values
    m1 = 100
    m2 = 1
    l = 2
    mu = 0.6
    g = 9.81
    # Integration region
    t_0 = 0
    t_end = 10

    args = [m1, m2, l, mu, g]
    initial_conds = [x0, phi0, v0, w0]

    sol = solve_ivp(diff_func, (t_0, t_end), initial_conds, method='BDF', args=([args]), dense_output=True)

    plot_solution_panels(sol)


if __name__ == '__main__':
    main()
