#!/usr/bin/env python3
import numpy as np


def diff_func_L1(t, vec, args):
    """
        Defines the differential equations for the pulley system using the Lagrangian approach.
        Method 1
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

    d_dt = [
        v,
        w,
        (m2 * x * (w ** 2) - g * (m1 + m2 * np.sin(phi))) / (m1 + m2) - +rhs,
        (-g * x * np.cos(phi) - 2 * x * v * w) / (x ** 2)  +rhs

    ]

    return d_dt


def diff_func_L2(t, vec, args):
    """
        Defines the differential equations for the pulley system using a Lagrangian approach.
        Method 2
        Arguments:
            vec :  vector of the state variables:
                      w = [x,phi,v,w]

            t :  time

            args : vector of the constant values:
                      p = [m1, m2, l, mu, g]

    """

    x, phi, v, w = vec
    m1, m2, l, mu, g = args
    try:
        rhs = (m1 * g) / np.exp(mu * -(phi + (np.pi / 2)))
    except:
        rhs = 10*20
    #rhs = (m1 * g) / np.exp(mu * -(phi + (np.pi / 2)))
    #print(rhs)

    l = l-(phi/2*np.pi)*0.2
    delta_x = l - x

    a = -m1 - m2
    b = -delta_x * m2 * np.sin(phi)
    c = m1 * g - m2 * v * w * np.sin(phi) - m2 * (w ** 2) * delta_x - m2 * g * np.sin(phi) + m2 * w * v * np.sin(
        phi) - m2 * (w ** 2) * delta_x * np.cos(phi) + rhs

    d = -m2 * delta_x * np.sin(phi)
    e = -m2 * (delta_x ** 2)
    f = m2 * g * delta_x * np.cos(phi) + m2 * (v ** 2) * np.sin(phi) + 2 * m2 * v * w * delta_x + rhs

    v_dot = (f * b - c * e) / (a * e - d * b)
    w_dot = (f * a - c * d) / (b * d - e * a)


    #print("w = {}, friction = {}".format(w, rhs))

    d_dt = [
        v,
        w,
        v_dot,
        w_dot
    ]

    return d_dt

