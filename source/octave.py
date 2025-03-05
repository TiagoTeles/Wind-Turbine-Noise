"""
Author:   T. Moreira da Fonte Fonseca Teles
Email:    tmoreiradafont@tudelft.nl
Date:     2025-03-05
License:  GNU GPL 3.0

Calculate the 1/3 octave band.

Classes:
    None

Functions:
    octave

Exceptions:
    None
"""

import numpy as np


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
        max_index =  np.ceil(10 * np.log10(f_max / f_ref) - 0.5)

        # Determine center, lower and upper frequencies
        f_center = f_ref * np.pow(10, np.arange(min_index, max_index+1) / 10)
        f_lower = f_center / np.pow(10, 1/20)
        f_upper = f_center * np.pow(10, 1/20)

    else:

        # Determine smallest and largest index
        min_index = np.floor(3 * np.log2(f_min / f_ref) + 0.5)
        max_index =  np.ceil(3 * np.log2(f_max / f_ref) - 0.5)

        # Determine center, lower and upper frequencies
        f_center = f_ref * np.pow(2, np.arange(min_index, max_index+1) / 3)
        f_lower = f_center / np.pow(2, 1/6)
        f_upper = f_center * np.pow(2, 1/6)

    return f_center, f_lower, f_upper
