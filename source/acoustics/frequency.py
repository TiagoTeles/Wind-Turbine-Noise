"""
Author:   T. Moreira da Fonte Fonseca Teles
Email:    tmoreiradafont@tudelft.nl
Date:     2025-11-10
License:  GNU GPL 3.0

Determine the one-third octave frequency bands.

Classes:
    None

Functions:
    one_third_octave
    doppler_factor

Exceptions:
    None
"""

import numpy as np

from source.constants import F_REF
from source.settings import F_MIN, F_MAX


def one_third_octave(f_min, f_max, base_10):
    """
    Determine the center, lower, and upper frequencies of the one-third octave frequency bands.

    Parameters:
        f_min : float -- minimum frequency, [Hz]
        f_max : float -- maximum frequency, [Hz]
        base_10 : bool -- use base-10 formulation?

    Returns:
        f_center : np.ndarray -- center frequencies, [Hz]
        f_lower : np.ndarray -- lower frequencies, [Hz]
        f_upper : np.ndarray -- upper frequencies, [Hz]
    """

    if base_10:

        # Determine the smallest and largest index
        min_index = np.floor(10.0 * np.log10(f_min / F_REF))
        max_index = np.ceil(10.0 * np.log10(f_max / F_REF))

        # Determine the center, lower and upper frequencies
        f_center = F_REF * np.pow(10.0, np.arange(min_index, max_index + 1) / 10.0)
        f_lower = f_center / np.pow(10.0, 1.0 / 20.0)
        f_upper = f_center * np.pow(10.0, 1.0 / 20.0)

    else:

        # Determine the smallest and largest index
        min_index = np.floor(3.0 * np.log2(f_min / F_REF))
        max_index = np.ceil(3.0 * np.log2(f_max / F_REF))

        # Determine the center, lower and upper frequencies
        f_center = F_REF * np.pow(2.0, np.arange(min_index, max_index + 1) / 3.0)
        f_lower = f_center / np.pow(2.0, 1.0 / 6.0)
        f_upper = f_center * np.pow(2.0, 1.0 / 6.0)

    return f_center, f_lower, f_upper


def doppler_factor(x_s, x_o, v_s, v_o, c_0):
    """
    Determine the frequency shift due to the Doppler effect.

    Parameters:
        x_s : np.ndarray -- source position, [m]
        x_o : np.ndarray -- observer position, [m]
        v_s : np.ndarray -- source velocity, [m/s]
        v_o : np.ndarray -- observer velocity, [m/s]
        c_0 : float -- speed of sound, [m/s]

    Returns:
        doppler_factor : np.ndarray -- f_o / f_s, [-]
    """

    # Determine the source-observer unit vector
    r_so = (x_o - x_s) / np.linalg.norm(x_o - x_s, axis=0)

    # Determine the velocity component in the direction of r_so
    v_s_r_so = np.sum(v_s * r_so, axis=0)
    v_o_r_so = np.sum(v_o * r_so, axis=0)

    # Determine the doppler factor
    doppler_factor = (c_0 - v_o_r_so) / (c_0 - v_s_r_so)

    return doppler_factor

if __name__ == "__main__":

    # Show the base-2 one-third octave frequency bands
    f_c, f_l, f_u = one_third_octave(F_MIN, F_MAX, False)
    print(f"Base-2 centre frequency:\n {f_c} [Hz] \n")

    # Show the base-10 one-third octave frequency bands
    f_c, f_l, f_u = one_third_octave(F_MIN, F_MAX, True)
    print(f"Base-10 centre frequency:\n {f_c} [Hz] \n")
