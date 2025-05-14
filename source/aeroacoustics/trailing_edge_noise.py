"""
Author:   T. Moreira da Fonte Fonseca Teles
Email:    tmoreiradafont@tudelft.nl
Date:     2025-05-12
License:  GNU GPL 3.0

Calculate the TBLTE noise spectra.

Classes:
    None

Functions:
    wall_pressure_spectrum
    spanwise_correlation_length
    E
    te_noise

Exceptions:
    None
"""

import numpy as np
import scipy as sp


def wall_pressure_spectrum(K_1, U, delta_star, rho_0):
    """
    Determine the wall pressure spectrum.

    # Arguments:
    #     K_1 : np.array -- streamwise wavenumber, [1/m]
    #     U : np.array -- velocity, [m/s]
    #     delta_star : np.array -- boundary layer displacement thickness, [m]
    #     rho_0 : float -- air density, [kg/m^3]

    # Returns:
    #     Phi_pp : np.array -- wall pressure spectrum, [Pa^2/Hz]
    """

    # Determine omega_tilde
    omega_tilde = K_1 * delta_star

    # Determine the spectrum function
    F = (33.28 * omega_tilde) / (1 - 5.489 * omega_tilde + 36.74 * np.square(omega_tilde) \
      + 0.1505 * np.pow(omega_tilde, 5))

    # Determine the wall pressure spectrum
    Phi_pp = np.square(0.5 * rho_0 * np.square(U)) * (delta_star / U) * 2E-5 * F

    return Phi_pp


def spanwise_correlation_length(omega, K_1, K_2, U_c, delta_star, b_c):
    """
    Determine the spanwise correlation length.

    Arguments:
        omega : np.array -- angular frequency, [rad/s]
        K_1 : float -- streamwise wavenumber, [1/m]
        K_2 : float -- spanwise wavenumber, [1/m]
        U_c : np.array -- convection velocity, [m/s]
        delta_star : np.array -- boundary layer displacement thickness, [m]
        b_c : float -- correlation coefficient, [-]

    Returns:
        l_y : np.array -- spanwise correlation length, [m]
    """

    # Determine omega_tilde
    omega_tilde = K_1 * delta_star

    # Determine the correlation length
    l_y = (omega / (b_c * U_c)) / (np.square(K_2) + np.square(omega) / np.square(b_c * U_c))

    # Apply the correction proposed by Casalino et al. (2019)
    l_y *= (1 - np.exp(-np.square(omega_tilde) / 0.09))

    return l_y


def E(x):
    """
    Determine the combination of Fresnel integrals.

    Arguments:
        x : np.array -- function argument, [-]

    Returns:
        E : np.array -- combination of Fresnel integrals, [-]
    """

    S_2, C_2 = sp.special.fresnel(np.sqrt(2 * x / np.pi))

    return C_2 - 1j * S_2


def te_noise(f, b, c, x, y, z, U, delta_star, c_0, rho_0, p_ref, alpha_c, b_c, base_10):
    """
    Determine the trailing-edge noise spectra.

    Arguments:
        f : np.array -- frequency, [Hz]
        b : np.array -- span, [m]
        c : np.array -- chord, [m]
        x : np.array -- x-coordinate, [m]
        y : np.array -- y-coordinate, [m]
        z : np.array -- z-coordinate, [m]
        U : np.array -- velocity, [m/s]
        delta_star : np.array -- boundary layer displacement thickness, [m]
        c_0 : float -- speed of sound, [m/s]
        rho_0 : float -- air density, [kg/m^3]
        p_ref : float -- reference pressure, [Pa]
        alpha_c : np.array -- speed ratio, [-]
        b_c : np.array -- correlation coefficient, [-]
        base_10 : bool -- Use base 10?
    """

    # Determine the Mach number
    M = U / c_0

    # Determine the Prandtl-Glauert factor
    beta = np.sqrt(1 - np.square(M))

    # Determine the convection velocity
    U_c = U / alpha_c

    # Determine distance corrected for convection effects
    S_0 = np.sqrt(np.square(x) + np.square(beta) * (np.square(y) + np.square(z)))

    # Determine the angular frequency
    omega = 2 * np.pi * f

    # Determine the acoustic wavenumber
    k = omega / c_0

    # Determine the convective wave number
    K = omega / U

    # Determine the aerodynamic wavenumbers
    K_1 = alpha_c * K
    K_2 = k * y / S_0

    # Determine the non-dimensional wavenumbers
    K_line = K * (c/2)
    K_1_line = K_1 * (c/2)
    K_2_line = K_2 * (c/2)

    # Determine the frequency parameters
    mu_line = K_line * M / np.square(beta)
    kappa_line = np.sqrt(np.square(mu_line) - np.square(K_2_line) / np.square(beta))

    # Determine the correction factor
    epsilon = np.pow(1 + 1 / (4 * mu_line), -0.5)

    # Determine the first-order main scattering term. Neglect the term np.exp(-2 * 1j * C)
    B = K_1_line + M * mu_line + kappa_line
    C = K_1_line - mu_line * (x / S_0 - M)

    f_1 = - np.exp(2 * 1j * C) / (1j * C) * ((1 + 1j) * np.exp(-2 * 1j * C) \
        * np.sqrt(B / (B - C)) * E(2 * (B - C)) - (1 + 1j) * E(2 * B) + 1)

    # Determine the second-order back-scattering correction
    D = kappa_line - mu_line * x / S_0

    G = (1 + epsilon) * np.exp(1j * (2 * kappa_line + D)) * np.sin(D - 2 * kappa_line) \
      / (D - 2 * kappa_line) + (1 - epsilon) * np.exp(1j * (-2 * kappa_line + D)) \
      * np.sin(D + 2 * kappa_line) / (D + 2 * kappa_line) + (1 + epsilon) * (1 - 1j) \
      / (2 * (D - 2 * kappa_line)) * np.exp(4 * 1j * kappa_line) * E(4 * kappa_line) \
      - (1 - epsilon) * (1 + 1j) / (2 * (D + 2 * kappa_line)) * np.exp(-4 * 1j * kappa_line) \
      * E(4 * kappa_line) + np.exp(2 * 1j * D) / 2 * np.sqrt(2 * kappa_line / D) * E(2 * D) \
      * ((1 + 1j) * (1 - epsilon) / (D + 2 * kappa_line) - (1 - 1j) * (1 + epsilon) \
      / (D - 2 * kappa_line))

    Theta = np.sqrt((K_1_line + mu_line * M + kappa_line) / (K_line + mu_line * M + kappa_line))

    H = (1 + 1j) * np.exp(-4 * 1j * kappa_line) * (1 - np.square(Theta)) \
      / (2 * np.sqrt(np.pi) * (alpha_c - 1) * K_line * np.sqrt(B))

    f_2 = np.exp(4 * 1j * kappa_line) * (1 - (1 + 1j) * E(4 * kappa_line))

    f_2.real = f_2.real
    f_2.imag = f_2.imag * epsilon

    f_2 = H * (f_2 - np.exp(2 * 1j * D) + 1j * (D + K_line + M * mu_line - kappa_line) * G)

    # Determine the radiation integral
    I = f_1 + f_2

    # Determine the wall pressure spectrum using Schlinker's model
    Phi_pp = wall_pressure_spectrum(K_1, U, delta_star, rho_0)

    # Determine the spanwise correlation length using Corcos' model
    l_y = spanwise_correlation_length(omega, K_1, K_2, U_c, delta_star, b_c)

    # Determine the streamwise-integrated wavenumber spectral density of wall-pressure ﬂuctuations
    Pi_0 = (1/np.pi) * Phi_pp * l_y

    # Determine the double-sided farfield acoustic PSD. Only valid if b >> c
    S_pp = np.square(omega * z * c / (4 * np.pi * c_0 * np.square(S_0))) * 2 * np.pi * b \
         * np.square(np.abs(I)) * Pi_0

    # Consider only one side of the airfoil
    # TODO: Check if this is correct
    # S_pp /= 4

    # Determine the bandwidth
    if base_10:
        f_lower = f / np.pow(10, 1/20)
        f_upper = f * np.pow(10, 1/20)
        delta_f = f_upper - f_lower

    else:
        f_lower = f / np.pow(2, 1/6)
        f_upper = f * np.pow(2, 1/6)
        delta_f = f_upper - f_lower

    # Determine the angular bandwidth
    delta_omega = 2 * np.pi * delta_f

    # Determine the 1/3 octave band SPL
    SPL = 10 * np.log10(2 * delta_omega * S_pp / np.square(p_ref))

    return SPL
