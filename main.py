#!/usr/bin/env python3
import numpy as np
from scipy.integrate import *
from plot import *


def diff_func(t, vec):
    """
        Defines the differential equations for the pulley system.

        Arguments:
            vec :  vector of the state variables:
                      w = [x,phi,v,w]
            t :  time

    """
    m2 = 2
    p = [12, m2, 2, 0.2, 9.81]
    x, phi, v, w = vec
    m1, m2, l, mu, g = p
    delta_x = l - x
    rhs = (m1 * g) / np.exp(mu * ((phi+(np.pi / 2))))
    # a = -m1
    # b = -m2
    # c = -m2 * delta_x * np.sin(phi)
    # d = m1 * g - m2 * v * w * np.sin(phi) - m2 * (w ** 2) * delta_x - m2 * g * np.sin(phi) + m2 * w * v * np.sin(
    #     phi) - m2 * (w ** 2) * delta_x * np.cos(phi) - rhs
    # e = -m2 * delta_x * np.sin(phi)
    # f = -m2 * (delta_x)**2
    # h = m2 * g * delta_x * np.cos(phi) + m2 * (v ** 2) * np.sin(phi) + 2 * m2 * v * w * delta_x - rhs

    a = -m1-m2
    b = -delta_x*m2*np.sin(phi)
    c = m1*g-m2*v*w*np.sin(phi)-m2*(w**2)*delta_x - m2*g*np.sin(phi)+m2*w*v*np.sin(phi) - m2*(w**2)*x*np.cos(phi) - rhs

    d = -m2*delta_x*np.sin(phi)
    e = -m2*(delta_x**2)
    f = m2*g*delta_x*np.cos(phi)+m2*(v**2)*np.sin(phi)+2*m2*v*w*delta_x-rhs
    d_dt = [
        v,
        w,
        (f*b-c*e)/(a*e-d*b),
        (f*a-c*d)/(b*d-e*a)
        # (d * f - g * c) / (f * (a + b) + c * e),
        # (h * (a + b) - d * e) / (c * e - f * (a + b))
    ]

    return d_dt


def main():
    x0 = 0
    v0 = 0
    phi0 = (np.pi) / 2
    w0 = 0

    initial_conds = [x0, phi0, v0, w0]
    sol = solve_ivp(diff_func, (0, 10), initial_conds, method='BDF', dense_output=True)

    plot_solution_panels(sol)


if __name__ == '__main__':
    main()
