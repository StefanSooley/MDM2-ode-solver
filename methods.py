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
    m1, m2, l, mu, g, r = args
    delta_x = l - x - phi*r

    rhs = np.exp(mu*phi)

    v_dot = (-rhs+m1*g+0.5*w*w*m2-m2*np.sin(phi))/(m1+m2)
    w_dot = (rhs-0.5*m2*w*w*r-r*m2*np.sin(phi)-phi*r*m2*np.cos(phi))/(m2*r)
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
    m1, m2, l, mu, g, r = args

    v_dot = (np.exp(mu * phi) * ((l-x) * w * w + g * np.sin(phi)) - ((m1 * g) / m2)) / (np.exp(mu * phi) + (m1 / m2))

    w_dot = (g * np.cos(phi) - 2 * v * w) / (l-x)
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
    m1, m2, l, mu, g, r = args

    if v < 0:
        # Stage 2, the acceleration and velocity will stay at 0.
        v = 0
        v_dot = 0

        # Since the rope has 0 thickness, m2 will continue to accelerate until it has Inf w which cannot be
        # calculated using floating point numbers. To fix this, the rotations will stop when the rope is 0.1m on the
        # LHS.

        # (l - x - phi*r) is the approximation for the length of the rope on the LHS

        if (l - x - phi * r) <= 0.01:

            # w must return to 0 instantaneously when m2 stops, so w_dot is set to -10^10000 to approximate an infinite
            # deceleration until w = 0.

            if w <= 0:
                w_dot = 0
                w = 0
            else:
                w_dot = -(10 * 10000)

        else:
            # w_dot equation while the length of the LHS is longer than 2m
            w_dot = (g * np.sin(phi) + r * w * w) / (l - x - phi * r)
    else:
        # Stage 1
        v_dot = (m1 * g + m2 * np.exp(mu * phi) * (g * np.cos(phi) - (l - (phi * r) - x) * (w ** 2))) / (
                m1 + m2 * np.exp(mu * phi))
        w_dot = (g * np.sin(phi) + 2 * v * w + r * w * w) / (l - phi * r - x)

    d_dt = [
        v,
        w,
        v_dot,
        w_dot
    ]

    return d_dt
