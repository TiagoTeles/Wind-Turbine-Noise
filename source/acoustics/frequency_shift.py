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


def frequency_shift(f_s, r, v_s, c_0):
    """
    Determine the frequency shift due to the Doppler effect.

    Parameters:
        f_s : np.array -- source frequency, [Hz]
        r : np.array -- source-observer vector, [m]
        v_s : np.array -- source velocity vector, [m/s]
        c_0 : float -- speed of sound, [m/s]

    Returns:
        f_o : np.array -- observer frequency, [Hz]
    """

    # Determine the dot product of r and v_s
    dot_product = r[0, :] * v_s[0, :] + r[1, :] * v_s[1, :] + r[2, :] * v_s[2, :]

    # Determine the norm of r and v_s
    r_norm = np.linalg.norm(r, axis=0)
    v_s_norm = np.linalg.norm(v_s, axis=0)

    # Determine the angle between r and v_s
    theta = np.arccos(dot_product / (r_norm * v_s_norm))

    # Determine the observer frequency
    f_o = f_s * c_0 / (c_0 - v_s_norm * np.cos(theta))

    return f_o
