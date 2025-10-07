"""
Author:   T. Moreira da Fonte Fonseca Teles
Email:    tmoreiradafont@tudelft.nl
Date:     2025-09-30
License:  GNU GPL 3.0

Determine the TBLTE noise spectra.

Classes:
    None

Functions:
    farfield_acoustic_psd
    E
    trailing_edge_noise

Exceptions:
    None
"""

import numpy as np
import scipy as sp

from boundary_layer.turbulence import spanwise_correlation_length, wall_pressure_spectrum


def farfield_acoustic_psd(f, b, c, U, Phi_pp, l_y, x, y, z, c_0, alpha_c):
    """
    Determine the far-field acoustic PSD.

    Parameters:
        f : np.ndarray -- frequency, [Hz]
        b : np.ndarray -- span, [m]
        c : np.ndarray -- chord, [m]
        U : np.ndarray -- velocity, [m/s]
        Phi_pp : np.ndarray -- wall pressure PSD, [Pa^2/Hz]
        l_y : np.ndarray -- spanwise correlation length, [m]
        x : np.ndarray -- x coordinate, [m]
        y : np.ndarray -- y coordinate, [m]
        z : np.ndarray -- z coordinate, [m]
        c_0 : float -- speed of sound, [m/s]
        alpha_c : np.ndarray -- speed ratio, [-]

    Returns:
        S_pp : np.ndarray -- farfield acoustic PSD, [Pa^2/Hz]
    """

    # Check for low aspect ratios
    if np.any(b/c < 3.0):
        print("Low aspect ratio detected! AR < 3 [-].")

    # Determine the Mach number
    M = U / c_0

    # Determine the Prandtl-Glauert factor
    beta = np.sqrt(1.0 - np.square(M))

    # Determine the distance corrected for convection effects
    S_0 = np.sqrt(np.square(x) + np.square(beta) * (np.square(y) + np.square(z)))

    # Determine the angular frequency
    omega = 2 * np.pi * f

    # Determine the acoustic wavenumber
    k = omega / c_0

    # Determine the convective wavenumber
    K = omega / U

    # Determine the aerodynamic wavenumbers
    K_1 = alpha_c * K
    K_2 = k * y / S_0

    # Determine the non-dimensional wavenumbers
    K_line = K * (c / 2.0)
    K_1_line = K_1 * (c / 2.0)
    K_2_line = K_2 * (c / 2.0)

    # Determine the frequency parameters
    mu_line = K_line * M / np.square(beta)
    kappa_line = np.sqrt(np.square(mu_line) - np.square(K_2_line) / np.square(beta))

    # Determine the correction factor
    epsilon = np.pow(1.0 + 1.0 / (4.0 * mu_line), -1.0 / 2.0)

    # Determine the first-order main scattering term. Neglect the term np.exp(-2 * 1j * C)
    B = K_1_line + M * mu_line + kappa_line
    C = K_1_line - mu_line * (x / S_0 - M)

    f_1 = - np.exp(2.0 * 1j * C) / (1j * C) * ((1.0 + 1j) * np.exp(-2.0 * 1j * C) \
        * np.sqrt(B / (B - C)) * E(2.0 * (B - C)) - (1.0 + 1j) * E(2.0 * B) + 1.0)

    # Determine the second-order back-scattering correction
    D = kappa_line - mu_line * x / S_0

    G = (1.0 + epsilon) * np.exp(1j * (2.0 * kappa_line + D)) * np.sin(D - 2.0 * kappa_line) \
      / (D - 2.0 * kappa_line) + (1.0 - epsilon) * np.exp(1j * (-2.0 * kappa_line + D)) \
      * np.sin(D + 2.0 * kappa_line) / (D + 2.0 * kappa_line) + (1.0 + epsilon) * (1.0 - 1j) \
      / (2.0 * (D - 2.0 * kappa_line)) * np.exp(4.0 * 1j * kappa_line) * E(4.0 * kappa_line) \
      - (1.0 - epsilon) * (1.0 + 1j) / (2.0 * (D + 2.0 * kappa_line)) * np.exp(-4.0 * 1j * kappa_line) \
      * E(4.0 * kappa_line) + np.exp(2.0 * 1j * D) / 2.0 * np.sqrt(2.0 * kappa_line / D) * E(2.0 * D) \
      * ((1.0 + 1j) * (1.0 - epsilon) / (D + 2.0 * kappa_line) - (1.0 - 1j) * (1.0 + epsilon) \
      / (D - 2.0 * kappa_line))

    Theta = np.sqrt((K_1_line + mu_line * M + kappa_line) / (K_line + mu_line * M + kappa_line))

    H = ((1.0 + 1j) * np.exp(-4.0 * 1j * kappa_line) * (1.0 - np.square(Theta))) \
      / (2.0 * np.sqrt(np.pi) * (alpha_c - 1.0) * K_line * np.sqrt(B))

    f_2 = np.exp(4.0 * 1j * kappa_line) * (1.0 - (1.0 + 1j) * E(4.0 * kappa_line))

    f_2.real = f_2.real
    f_2.imag = f_2.imag * epsilon

    f_2 = H * (f_2 - np.exp(2.0 * 1j * D) + 1j * (D + K_line + M * mu_line - kappa_line) * G)

    # Determine the radiation integral
    I = f_1 + f_2

    # Determine the stream-wise-integrated wavenumber spectral density of wall-pressure ﬂuctuations
    Pi_0 = (1.0 / np.pi) * Phi_pp * l_y

    # Determine the double-sided farfield acoustic PSD
    S_pp = np.square((omega * z * c) / (4.0 * np.pi * c_0 * np.square(S_0))) \
         * 2.0 * np.pi * b * np.square(np.abs(I)) * Pi_0

    return S_pp


def E(x):
    """
    Determine the Fresnel integral.

    Parameters:
        x : np.ndarray -- function argument, [-]

    Returns:
        E : np.ndarray -- Fresnel integral, [-]
    """

    S_2, C_2 = sp.special.fresnel(np.sqrt(2 * x / np.pi))

    return C_2 - 1j * S_2


def trailing_edge_noise(f, b, c, U, delta_star, x, y, z, c_0, rho_0, b_c, alpha_c, p_ref):
    """
    Determine the TBLTE noise SPL.

    Parameters:
        f : np.ndarray -- frequency, [Hz]
        b : np.ndarray -- span, [m]
        c : np.ndarray -- chord, [m]
        U : np.ndarray -- velocity, [m/s]
        delta_star : np.ndarray -- boundary layer displacement thickness, [m]
        x : np.ndarray -- x coordinate, [m]
        y : np.ndarray -- y coordinate, [m]
        z : np.ndarray -- z coordinate, [m]
        c_0 : float -- speed of sound, [m/s]
        rho_0 : float -- air density, [kg/m^3]
        b_c : np.ndarray -- correlation coefficient, [-]
        alpha_c : np.ndarray -- speed ratio, [-]
        p_ref : float -- reference pressure, [Pa]

    Returns:
        spl : np.ndarray -- TBLTE noise SPL, [dB]
    """

    # Determine the wall pressure spectrum using Schlinker's model
    Phi_pp = wall_pressure_spectrum(f, U, delta_star, rho_0)

    # Determine the angular frequency
    omega = 2.0 * np.pi * f

    # Determine the acoustic wavenumber
    k = omega / c_0

    # Determine the Mach number
    M = U / c_0

    # Determine the Prandtl-Glauert factor
    beta = np.sqrt(1.0 - np.square(M))

    # Determine the distance corrected for convection effects
    S_0 = np.sqrt(np.square(x) + np.square(beta) * (np.square(y) + np.square(z)))

    # Determine the spanwise wavenumber
    K_2 = k * y / S_0

    # Determine the spanwise correlation length using Corcos' model
    l_y = spanwise_correlation_length(f, U, delta_star, K_2, b_c, alpha_c)

    # Determine the far-field acoustic PSD
    S_pp = farfield_acoustic_psd(f, b, c, U, Phi_pp, l_y, x, y, z, c_0, alpha_c)

    # Determine the bandwidth
    delta_omega = 2.0 * np.pi * (0.231 * f)

    # Determine the 1/3 octave band SPL
    spl = 10.0 * np.log10(2.0 * S_pp * delta_omega / np.square(p_ref))

    return spl
