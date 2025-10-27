"""
Author:   T. Moreira da Fonte Fonseca Teles
Email:    tmoreiradafont@tudelft.nl
Date:     2025-10-06
License:  GNU GPL 3.0

Determine the octave frequency bands.

Classes:
    None

Functions:
    one_third_octave
    doppler_effect

Exceptions:
    None
"""

import numpy as np

from source.settings import F_MIN, F_MAX, F_REF


def one_third_octave(f_min, f_max, f_ref, base_10):
    """
    Determine the center, lower, and upper frequencies of the 1/3 octave frequency bands.

    Parameters:
        f_min : float -- minimum frequency, [Hz]
        f_max : float -- maximum frequency, [Hz]
        f_ref : float -- reference frequency, [Hz]
        base_10 : bool -- use base 10?

    Returns:
        f_center : np.ndarray -- center frequencies, [Hz]
        f_lower : np.ndarray -- lower frequencies, [Hz]
        f_upper : np.ndarray -- upper frequencies, [Hz]
    """

    if base_10:

        # Determine the smallest and largest index
        min_index = np.floor(10.0 * np.log10(f_min / f_ref))
        max_index = np.ceil(10.0 * np.log10(f_max / f_ref))

        # Determine the center, lower and upper frequencies
        f_center = f_ref * np.pow(10.0, np.arange(min_index, max_index + 1) / 10.0)
        f_lower = f_center / np.pow(10.0, 1.0 / 20.0)
        f_upper = f_center * np.pow(10.0, 1.0 / 20.0)

    else:

        # Determine the smallest and largest index
        min_index = np.floor(3.0 * np.log2(f_min / f_ref))
        max_index = np.ceil(3.0 * np.log2(f_max / f_ref))

        # Determine the center, lower and upper frequencies
        f_center = f_ref * np.pow(2.0, np.arange(min_index, max_index + 1) / 3.0)
        f_lower = f_center / np.pow(2.0, 1.0 / 6.0)
        f_upper = f_center * np.pow(2.0, 1.0 / 6.0)

    return f_center, f_lower, f_upper


def doppler_effect(f_s, x_s, x_o, v_s, v_o, c_0):
    """
    Determine the frequency shift due to the Doppler effect.

    Parameters:
        f_s : np.ndarray -- source frequency, [Hz]
        x_s : np.ndarray -- source position, [m]
        x_o : np.ndarray -- observer position, [m]
        v_s : np.ndarray -- source velocity, [m/s]
        v_o : np.ndarray -- observer velocity, [m/s]
        c_0 : float -- speed of sound, [m/s]

    Returns:
        f_o : np.ndarray -- observer frequency, [Hz]
    """

    # Determine the source-observer unit vector
    r_so = (x_o - x_s) / np.linalg.norm(x_o - x_s, axis=0)

    # Determine the velocity component in the direction of r_so
    v_s_r_so = np.sum(v_s * r_so, axis=0)
    v_o_r_so = np.sum(v_o * r_so, axis=0)

    # Determine the observer frequency
    f_o = f_s * (c_0 - v_o_r_so) / (c_0 - v_s_r_so)

    return f_o

if __name__ == "__main__":

    # Show the base-2 one-third octave frequency bands
    f_c, f_l, f_u = one_third_octave(F_MIN, F_MAX, F_REF, False)
    print("Base-2 centre frequency:\n" + str(f_c) + " [Hz] \n")

    # Show the base-10 one-third octave frequency bands
    f_c, f_l, f_u = one_third_octave(F_MIN, F_MAX, F_REF, True)
    print("Base-10 centre frequency:\n" + str(f_c) + " [Hz] \n")
