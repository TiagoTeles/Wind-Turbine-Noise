"""
Author:   T. Moreira da Fonte Fonseca Teles
Email:    tmoreiradafont@tudelft.nl
Date:     2025-05-05
License:  GNU GPL 3.0

Miscellaneous functions.

Classes:
    None

Functions:
    octave
    turbulence_intensity
    turbulence_length_scale
    surface_roughness_length

Exceptions:
    None
"""

import numpy as np


def octave(f_min, f_max, f_ref, base_10):
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


def turbulence_intensity(z, z_0):
    """
    Determine the turbulence intensity.

    Parameters:
        z : np.array -- height above the terrain, [m]
        z_0 :  np.array -- surface roughness length, [m]

    Returns:
        I :  np.array -- turbulence intensity, [-]
    """

    gamma = 0.24 + 0.096 * np.log10(z_0) + 0.016 * np.square(np.log10(z_0))
    I = gamma * np.log10(30/z_0) / np.log10(z/z_0)

    return I


def turbulence_length_scale(z, z_0):
    """
    Determine the turbulence length scale.

    Parameters:
        z : np.array -- height above the terrain, [m]
        z_0 :  np.array -- surface roughness length, [m]

    Returns:
        L :  np.array -- turbulence length scale, [m]
    """

    L = 25 * np.pow(z, 0.35) * np.pow(z_0, -0.063)

    return L


def surface_roughness_length(z_ref, U_ref, g, kappa, nu):
    """
    Determine the surface roughness length.

    Parameters:
        z_ref : np.array -- reference height, [m]
        U_ref : np.array -- reference velocity, [m/s]
        g : np.array -- gravitational acceleration, [m/s^2]
        kappa : np.array -- Von Karman constant, [-]
        nu : np.array -- kinematic viscosity, [m^2/s]

    Returns:
        z_0 : np.array -- surface roughness length, [m]
    """

    # Determine the R and A
    R = z_ref / (0.11 * nu) * (kappa * U_ref)
    A = 0.018 / (g * z_ref) * np.square(kappa * U_ref)

    # Determine b_n
    b_n_nu = -1.47 + 0.93 * np.log(R)
    b_n_alpha = 2.65 - 1.44 * np.log(A) - 0.015 * np.square(np.log(A))
    b_n = np.pow(np.pow(b_n_nu, -12) + np.pow(b_n_alpha, -12), -1/12)

    # Determine the roughness length
    z_0 = z_ref / (np.exp(b_n) - 1)

    return z_0
