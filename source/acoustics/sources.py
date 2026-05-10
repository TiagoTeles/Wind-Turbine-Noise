"""
Author:   T. Moreira da Fonte Fonseca Teles
Email:    tmoreiradafont@tudelft.nl
Date:     2026-05-10
License:  GNU GPL 3.0

Determine acoustic multipole contributions.

Classes:
    None

Functions:
    monopole
    dipole

Exceptions:
    None
"""

import numpy as np


def monopole(f, x_s, x_o, c_0):
    """
    Determine the complex pressure of a monopole.

    Parameters:
        f : np.ndarray -- frequency, [Hz]
        x_s : np.ndarray -- source position, [m]
        x_o : np.ndarray -- observer position, [m]
        c_0 : float -- speed of sound, [m/s]

    Returns:
        p : np.ndarray -- pressure, [Pa]
    """

    # Determine the wavenumber
    k = 2.0 * np.pi * f / c_0

    # Determine the distance from the source to the observer
    r_so = np.linalg.norm(x_o - x_s, axis=0)

    # Determine the complex pressure
    p = np.exp(1.0j * k * r_so) / (4.0 * np.pi * r_so)

    return p


def dipole(f, x_s, x_o, n_s, d, c_0):
    """
    Determine the complex pressure of a dipole.

    Parameters:
        f : np.ndarray -- frequency, [Hz]
        x_s : np.ndarray -- source position, [m]
        x_o : np.ndarray -- observer position, [m]
        n_s : np.ndarray -- dipole axis, [-]
        d : np.ndarray -- monopole distance, [m]
        c_0 : float -- speed of sound, [m/s]

    Returns:
        p : np.ndarray -- pressure, [Pa]
        x_s_0 : np.ndarray -- source 0 position, [m]
        x_s_1 : np.ndarray -- source 1 position, [m]
    """

    # Ensure the axis direction is normalised
    n_s[0, :, :] /= np.linalg.norm(n_s, axis=0)
    n_s[1, :, :] /= np.linalg.norm(n_s, axis=0)
    n_s[2, :, :] /= np.linalg.norm(n_s, axis=0)

    # Determine the monopole positions
    x_s_0 = x_s + n_s * (d / 2.0)
    x_s_1 = x_s - n_s * (d / 2.0)

    # Determine the monopole pressures
    p_0  = monopole(f, x_s_0, x_o, c_0)
    p_1  = monopole(f, x_s_1, x_o, c_0)

    # Subtract the monopole pressures
    p = p_0 - p_1

    return p, x_s_0, x_s_1