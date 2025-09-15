"""
Author:   T. Moreira da Fonte Fonseca Teles
Email:    tmoreiradafont@tudelft.nl
Date:     2025-09-05
License:  GNU GPL 3.0

Determine the frequency shift due to the Doppler effect.

Classes:
    None

Functions:
    frequency_shift

Exceptions:
    None
"""

import numpy as np


def frequency_shift(f_s, x_s, x_o, v_s, v_o, c_0):
    """
    Determine the frequency shift due to the Doppler effect.

    Parameters:
        f_s : np.array -- source frequency, [Hz]
        x_s : np.array -- source position, [m]
        x_o : np.array -- observer position, [m]
        v_s : np.array -- source velocity, [m/s]
        v_o : np.array -- observer velocity, [m/s]
        c_0 : float -- speed of sound, [m/s]

    Returns:
        f_o : np.array -- observer frequency, [Hz]
    """

    # Determine the source-observer unit vector
    r_so = (x_o - x_s) / np.linalg.norm(x_o - x_s, axis=0)

    # Determine the velocity in the direction of r_so
    v_s_proj = np.sum(v_s * r_so, axis=0)
    v_o_proj = np.sum(v_o * r_so, axis=0)

    # Determine the observer frequency
    f_o = f_s * (c_0 - v_o_proj) / (c_0 - v_s_proj)

    return f_o
