"""
Author:   T. Moreira da Fonte Fonseca Teles
Email:    tmoreiradafont@tudelft.nl
Date:     2025-02-21
License:  GNU GPL 3.0

Calculate the 1/3 octave band.

Classes:
    None

Functions:
    None

Exceptions:
    None
"""

import numpy as np


# Settings
BASE_2 = False  # Base 2 or 10?
F_MIN = 20      # Minimum frequency, [Hz]
F_MAX = 20000   # Maximum frequency, [Hz]
F_REF = 1000    # Reference frequency, [Hz]

if BASE_2:
    # Determine smallest and largest index
    min_index = np.floor(3 * np.log2(F_MIN / F_REF) + 0.5)
    max_index =  np.ceil(3 * np.log2(F_MAX / F_REF) - 0.5)

    # Determine center, lower and upper frequencies
    f_center = F_REF * np.pow(2, np.arange(min_index, max_index+1) / 3)
    f_lower = f_center / np.pow(2, 1/6)
    f_upper = f_center * np.pow(2, 1/6)

else:
    # Determine smallest and largest index
    min_index = np.floor(10 * np.log10(F_MIN / F_REF) + 0.5)
    max_index =  np.ceil(10 * np.log10(F_MAX / F_REF) - 0.5)

    # Determine center, lower and upper frequencies
    f_center = F_REF * np.pow(10, np.arange(min_index, max_index+1) / 10)
    f_lower = f_center / np.pow(10, 1/20)
    f_upper = f_center * np.pow(10, 1/20)
