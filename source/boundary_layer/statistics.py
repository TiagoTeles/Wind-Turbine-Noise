"""
Author:   T. Moreira da Fonte Fonseca Teles
Email:    tmoreiradafont@tudelft.nl
Date:     2025-11-11
License:  GNU GPL 3.0

Determine the boundary layer statistics.

Classes:
    None

Functions:
    wall_pressure_spectrum
    spanwise_correlation_length

Exceptions:
    None
"""

import numpy as np

from source.constants import ALPHA_C, B_C


def wall_pressure_spectrum(f, U, delta_star, rho):
    """
    Determine the wall pressure spectrum.

    Parameters:
        f : np.ndarray -- frequency, [Hz]
        U : np.ndarray -- velocity, [m/s]
        delta_star : np.ndarray -- displacement thickness, [m]
        rho : float -- density, [kg/m^3]

    Returns:
        Phi_pp : np.ndarray -- wall pressure spectrum, [Pa^2/Hz]
    """

    # Determine the angular frequency
    omega = 2.0 * np.pi * f

    # Determine the convective wavenumber
    K = omega / U

    # Determine the non-dimensional wavenumber
    omega_tilde = K * delta_star

    # Determine the spectrum function
    F = (33.28 * omega_tilde) / (1.0 - 5.489 * omega_tilde + 36.74 * np.square(omega_tilde) \
      + 0.1505 * np.pow(omega_tilde, 5.0))

    # Determine the wall pressure spectrum
    Phi_pp = np.square(0.5 * rho * np.square(U)) * (delta_star / U) * 2.0E-5 * F

    return Phi_pp


def spanwise_correlation_length(f, U, delta_star, K_2):
    """
    Determine the spanwise correlation length.

    Parameters:
        f : np.ndarray -- frequency, [Hz]
        U : np.ndarray -- velocity, [m/s]
        delta_star : np.ndarray -- displacement thickness, [m]
        K_2 : np.ndarray -- spanwise aerodynamic wavenumber, [1/m]

    Returns:
        l_y : np.ndarray -- spanwise correlation length, [m]
    """

    # Determine the angular frequency
    omega = 2.0 * np.pi * f

    # Determine the convective velocity
    U_c = U / ALPHA_C

    # Determine the spanwise correlation length
    l_y = (omega / (B_C * U_c)) / (np.square(K_2) + np.square(omega) / np.square(B_C * U_c))

    # Determine the convective wavenumber
    K = omega / U

    # Determine the non-dimensional wavenumber
    omega_tilde = K * delta_star

    # Apply the correction proposed by Casalino et al. (2022)
    l_y *= (1.0 - np.exp(-np.square(omega_tilde) / 0.09))

    return l_y
