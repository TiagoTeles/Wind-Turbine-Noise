"""
Author:   T. Moreira da Fonte Fonseca Teles
Email:    tmoreiradafont@tudelft.nl
Date:     2025-08-28
License:  GNU GPL 3.0

Determine the octave frequency bands.

Classes:
    None

Functions:
    one_third_octave

Exceptions:
    None
"""

import numpy as np


def one_third_octave(f_min, f_max, f_ref, base_10):
    """
    Determine the center, lower, and upper frequencies of the 1/3 octave frequency bands.

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
        # Determine the smallest and largest index
        min_index = np.floor(10 * np.log10(f_min / f_ref))
        max_index = np.ceil(10 * np.log10(f_max / f_ref))

        # Determine the center, lower and upper frequencies
        f_center = f_ref * np.pow(10, np.arange(min_index, max_index+1) / 10)
        f_lower = f_center / np.pow(10, 1/20)
        f_upper = f_center * np.pow(10, 1/20)

    else:
        # Determine the smallest and largest index
        min_index = np.floor(3 * np.log2(f_min / f_ref))
        max_index = np.ceil(3 * np.log2(f_max / f_ref))

        # Determine the center, lower and upper frequencies
        f_center = f_ref * np.pow(2, np.arange(min_index, max_index+1) / 3)
        f_lower = f_center / np.pow(2, 1/6)
        f_upper = f_center * np.pow(2, 1/6)

    return f_center, f_lower, f_upper

if __name__ == "__main__":

    # Show the base-2 one-third octave frequency bands
    f_c, f_l, f_u = one_third_octave(20, 20000, 1000, base_10=False)
    print("Base 2:\n" + str(f_c) + " [Hz] \n")

    # Show the base-10 one-third octave frequency bands
    f_c, f_l, f_u = one_third_octave(20, 20000, 1000, base_10=True)
    print("Base 10:\n" + str(f_c) + " [Hz] \n")
