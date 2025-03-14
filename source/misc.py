"""
Author:   T. Moreira da Fonte Fonseca Teles
Email:    tmoreiradafont@tudelft.nl
Date:     2025-03-14
License:  GNU GPL 3.0

Miscellaneous functions.

Classes:
    None

Functions:
    octave
    E

Exceptions:
    None
"""

import numpy as np
import scipy as sp


def octave(f_min, f_max, f_ref, base_10=True):
    """
    Determine the center, lower, and upper frequencies of the 1/3 octave band.

    Arguments:
        f_min : float -- Minimum frequency, [Hz]
        f_max : float -- Maximum frequency, [Hz]
        f_ref : float -- Reference frequency, [Hz]
        base_10 : bool -- Use base 10?

    Returns:
        f_center : np.array -- Center frequencies, [Hz]
        f_lower : np.array -- Lower frequencies, [Hz]
        f_upper : np.array -- Upper frequencies, [Hz]
    """

    if base_10:
        # Determine smallest and largest index
        min_index = np.floor(10 * np.log10(f_min / f_ref) + 0.5)
        max_index = np.ceil(10 * np.log10(f_max / f_ref) - 0.5)

        # Determine center, lower and upper frequencies
        f_center = f_ref * np.pow(10, np.arange(min_index, max_index+1) / 10)
        f_lower = f_center / np.pow(10, 1/20)
        f_upper = f_center * np.pow(10, 1/20)

    else:
        # Determine smallest and largest index
        min_index = np.floor(3 * np.log2(f_min / f_ref) + 0.5)
        max_index = np.ceil(3 * np.log2(f_max / f_ref) - 0.5)

        # Determine center, lower and upper frequencies
        f_center = f_ref * np.pow(2, np.arange(min_index, max_index+1) / 3)
        f_lower = f_center / np.pow(2, 1/6)
        f_upper = f_center * np.pow(2, 1/6)

    return f_center, f_lower, f_upper


def E(x):
    """
    Determine the combination of Fresnel integrals defined in Roger and Moreau (2005).

    Arguments:
        x : np.array -- function argument, [-]

    Returns:
        E : np.array -- combination of Fresnel integrals, [-]
    """

    s_2, c_2 = sp.special.fresnel(np.sqrt(2 * x / np.pi))

    return c_2 - 1j * s_2


def sears(x):
    """
    Sears function approximation.

    Arguments:
        x : np.array -- function argument, [-]

    Returns:
        S : np.array -- Sears function, [-]
    """

    return np.sqrt(1 / (2*np.pi*x + 1 / (1 + 2.4*x)))
