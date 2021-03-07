#!/usr/bin/env python3
import numpy as np


def diff_func_L1(t, vec, args):
    """
        Defines the differential equations for the pulley system using the Lagrangian approach.
        Overleaf set of equations.
        Arguments:
            vec :  vector of the state variables:
                      w = [x,phi,v,w]

            t :  time

            args : vector of the constant values:
                      p = [m1, m2, l, mu, g]

    """

    x, phi, v, w = vec
    m1, m2, l, mu, g = args

    rhs = (m1 * g) / np.exp(mu * -(phi + (np.pi / 2)))
    delta_x = l - x

    d_dt = [
        v,
        w,
        (m2 * x * (w ** 2) - g * (m1 + m2 * np.sin(phi))) / (m1 + m2) + rhs,
        (-g * x * np.cos(phi) - 2 * x * v * w) / (x ** 2) + rhs

    ]

    return d_dt


def diff_func_L2(t, vec, args):
    """
        Defines the differential equations for the pulley system using a Lagrangian approach.
        Original set of equations
        Arguments:
            vec :  vector of the state variables:
                      w = [x,phi,v,w]

            t :  time

            args : vector of the constant values:
                      p = [m1, m2, l, mu, g]

    """

    x, phi, v, w = vec
    m1, m2, l, mu, g = args

    rhs = (m1 * g) / np.exp(mu * -(phi + (np.pi / 2)))

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

    # print("w = {}, friction = {}".format(w, rhs))

    d_dt = [
        v,
        w,
        v_dot,
        w_dot
    ]

    return d_dt


def diff_func_L3(t, vec, args):
    """
        Defines the differential equations for the pulley system using a Lagrangian approach.
        2nd set of equations
        Arguments:
            vec :  vector of the state variables:
                      w = [x,phi,v,w]

            t :  time

            args : vector of the constant values:
                      p = [m1, m2, l, mu, g]

    """

    x, phi, v, w = vec
    m1, m2, l, mu, g = args

    delta_x = l - x

    rhs = (m1 * g) / np.exp(mu * -(phi + (np.pi / 2)))

    v_dot = (rhs - g * m1 * np.sin(phi) - 2 * m2 * (w ** 2) * (delta_x) - g * m2) / m2
    w_dot = (rhs + 4 * l * m1 * v * w + m2 * g * delta_x * np.sin(phi) - 4 * m1 * x * v * w) / (
            2 * (l ** 2) * m1 + 2 * m1 * (x ** 2) - 4 * l * m1 * x)

    d_dt = [
        v,
        w,
        v_dot,
        w_dot
    ]

    return d_dt


def diff_func_L4(t, vec, args):
    """
        Defines the differential equations for the pulley system using a Lagrangian approach.
        Uses the overleaf lagrangian equations
        Arguments:
            vec :  vector of the state variables:
                      w = [x,phi,v,w]

            t :  time

            args : vector of the constant values:
                      p = [m1, m2, l, mu, g]

    """

    x, phi, v, w = vec
    m1, m2, l, mu, g = args

    rhs = (m1 * g) / np.exp(mu * -(phi + (np.pi / 2)))
    # rhs = mu*phi
    v_dot = (m2 * x * w * w - g * (m1 + m2 * np.sin(phi)) - rhs) / (m1 + m2)
    w_dot = (-rhs - (g * x * m2 * np.cos(phi)) - 2 * m2 * w * x * v) / (m2 * x * x)

    d_dt = [
        v,
        w,
        v_dot,
        w_dot
    ]

    return d_dt


def diff_func_L5(t, vec, args):
    """
        Defines the differential equations for the pulley system using a Lagrangian approach.
        Uses the overleaf lagrangian equations
        Arguments:
            vec :  vector of the state variables:
                      w = [x,phi,v,w]

            t :  time

            args : vector of the constant values:
                      p = [m1, m2, l, mu, g]

    """

    x, phi, v, w = vec
    m1, m2, l, mu, g = args
    delta_x = l - x

    rhs = 1 / np.exp(mu * -(phi + (np.pi / 2)))

    v_dot = (rhs - m1 * g + m2 * g * np.sin(phi) + m2 * (delta_x) * w * w) / (m1 + m2)
    w_dot = ((rhs / m2 * g) + (delta_x) * np.cos(phi)) / (x * x - 4 * l * x)
    d_dt = [
        v,
        w,
        v_dot,
        w_dot
    ]

    return d_dt


def diff_func_N1(t, vec, args):
    """
        Defines the differential equations for the pulley system using a Newtonian approach.
        Uses the overleaf newtonian equations
        Arguments:
            vec :  vector of the state variables:
                      w = [x,phi,v,w]

            t :  time

            args : vector of the constant values:
                      p = [m1, m2, l, mu, g]

    """

    x, phi, v, w = vec
    m1, m2, l, mu, g = args

    v_dot = (np.exp(mu * phi) * (x * w * w + g * np.sin(phi)) - ((m1 * g) / m2)) / (np.exp(mu * phi) + (m1 / m2))
    print(v_dot)
    w_dot = (g * np.cos(phi) - 2 * v * w) / x
    d_dt = [
        v,
        w,
        v_dot,
        w_dot
    ]

    return d_dt


def diff_func_P(t, vec, args):
    """
        Defines the differential equations for the pulley system using a Newtonian approach.
        Uses the equations from the paper
        Arguments:
            vec :  vector of the state variables:
                      w = [x,phi,v,w]

            t :  time

            args : vector of the constant values:
                      p = [m1, m2, l, mu, g]

    """

    x, phi, v, w = vec
    m1, m2, l, mu, g = args
    r = 1  # Radius of the pulley in m

    if v < 0:
        v_dot = 0
        if l - x - phi * r <= 0.3:
            w,w_dot = 0,0

        else:
            w_dot = (-g * np.sin(phi) + r * w * w) / (l - x - phi * r)
    else:
        v_dot = (m1 * g + m2 * np.exp(mu * phi) * (g * np.cos(phi) - (l - (phi * r) - x) * (w ** 2))) / (
                m1 + m2 * np.exp(mu * phi))
        w_dot = (g * np.sin(phi) + 2 * v * w + x * w * w) / (l - phi * r - x)

    d_dt = [
        v,
        w,
        v_dot,
        w_dot
    ]
    print(d_dt)
    return d_dt
