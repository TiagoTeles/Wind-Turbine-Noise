"""
Author:   T. Moreira da Fonte Fonseca Teles
Email:    tmoreiradafont@tudelft.nl
Date:     2025-09-05
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

from source.settings import BASE_10, F_MAX, F_MIN, F_REF


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
        min_index = np.floor(10.0 * np.log10(f_min / f_ref))
        max_index = np.ceil(10.0 * np.log10(f_max / f_ref))

        # Determine the center, lower and upper frequencies
        f_center = f_ref * np.pow(10.0, np.arange(min_index, max_index + 1) / 10.0)
        f_lower = f_center / np.pow(10.0, 1.0/20.0)
        f_upper = f_center * np.pow(10.0, 1.0/20.0)

    else:

        # Determine the smallest and largest index
        min_index = np.floor(3.0 * np.log2(f_min / f_ref))
        max_index = np.ceil(3.0 * np.log2(f_max / f_ref))

        # Determine the center, lower and upper frequencies
        f_center = f_ref * np.pow(2.0, np.arange(min_index, max_index + 1) / 3.0)
        f_lower = f_center / np.pow(2.0, 1.0/6.0)
        f_upper = f_center * np.pow(2.0, 1.0/6.0)

    return f_center, f_lower, f_upper

if __name__ == "__main__":

    # Show the base-2 one-third octave frequency bands
    f_c, f_l, f_u = one_third_octave(F_MIN, F_MAX, F_REF, False)
    print("Base-2 centre frequency:\n" + str(f_c) + " [Hz] \n")

    # Show the base-10 one-third octave frequency bands
    f_c, f_l, f_u = one_third_octave(F_MIN, F_MAX, F_REF, True)
    print("Base-2 centre frequency:\n" + str(f_c) + " [Hz] \n")
