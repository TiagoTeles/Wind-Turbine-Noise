"""
Author:   T. Moreira da Fonte Fonseca Teles
Email:    tmoreiradafont@tudelft.nl
Date:     2025-03-24
License:  GNU GPL 3.0

Calculate the TBLTE noise spectra.

Classes:
    None

# Functions:
    delta_star
    wall_pressure
    correlation_length
    te_noise

Exceptions:
    None
"""

import os
import sys

import numpy as np

from misc import E
from settings import XFOIL_PATH
from xfoil import XFoil


def displacement_thickness(blade, pos, Re, M, alpha, cutoff, probe_top, probe_bot, it):
    """
    Determine the boundary layer displacement thickness.

    Arguments:
        blade : Blade -- blade object
        pos : np.array -- position along the blade, [-]
        Re : np.array -- Reynolds number, [-]
        M : np.array -- Mach number, [-]
        alpha : np.array -- angle of attack, [rad]
        cutoff : float -- radial cutoff for the airfoil, [-]
        probe_top : np.array -- probe location at the top surface, [-]
        probe_bot : np.array -- probe location at the bottom surface, [-]
        it : int -- number of iterations for XFoil

    Returns:
        delta_star_top : np.array -- top boundary layer displacement thickness, [m]
        delta_star_bot : np.array -- bottom boundary layer displacement thickness, [m]
    """

    # Remove the unused axis
    Re = np.squeeze(Re)
    M = np.squeeze(M)
    alpha = np.squeeze(alpha)

    # Initialise the numpy arrays
    delta_star_top = np.zeros(Re.shape)
    delta_star_bot = np.zeros(Re.shape)

    # Iterate through each panel
    for i in range(pos):

        airfoil = blade.interpolate(pos[i])

        # Use a cutoff to ensure that XFoil converges
        # For cutoff=0.4, SPL_cutoff = SPL_tip - 19.90 dB
        if pos[i] / blade.data["pos"].iloc[-1]  > cutoff:

            # Determine the path to the airfoil
            cwd = os.path.dirname(blade.path)
            path = os.path.relpath(airfoil.path, cwd)

            # Run XFoil
            xfoil = XFoil(XFOIL_PATH, cwd)
            top, bot = xfoil.run(path, Re[i], M[i], np.degrees(alpha[i]), it=it)

            # Determine the displacement thickness at the probe locations
            delta_star_top[i] = np.interp(probe_top, top["x/c"], top["delta_star"])
            delta_star_bot[i] = np.interp(probe_bot, bot["x/c"], bot["delta_star"])

        else:

            # Set the displacement thickness to zero
            delta_star_top[i] = 0
            delta_star_bot[i] = 0

    return delta_star_top, delta_star_bot


def wall_pressure(K_1, U, delta_star, rho_0):
    """
    Determine the wall pressure spectrum.

    Arguments:
        K_1 : float -- streamwise wavenumber, [1/m]
        U : np.array -- velocity, [m/s]
        delta_star : np.array -- boundary layer displacement thickness, [m]
        rho_0 : float -- air density, [kg/m^3]

    Returns:
        Phi_pp : np.array -- wall pressure spectrum, [Pa^2/Hz]
    """

    omega_tilde = K_1 * delta_star

    F = (33.28 * omega_tilde) / (1 - 5.489 * omega_tilde + 36.74 * np.square(omega_tilde) \
      + 0.1505 * np.pow(omega_tilde, 5))

    Phi_pp = np.square(0.5 * rho_0 * np.square(U)) * (delta_star / U) * 2E-5 * F

    return Phi_pp


def correlation_length(omega, K_1, K_2, U_c, delta_star, b_c):
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


def te_noise(f, blade, pos, b, c, x, y, z, U, Re, alpha, alpha_c, b_c, p_ref, c_0, rho_0, cutoff, \
             probe_top, probe_bot, it):
    """
    Determine the trailing-edge noise spectra.

    Arguments:
        f : np.array -- frequency, [Hz]
        blade : Blade -- blade object
        pos : np.array -- position along the blade, [-]
        b : np.array -- span, [m]
        c : np.array -- chord, [m]
        x : np.array -- x-coordinate, [m]
        y : np.array -- y-coordinate, [m]
        z : np.array -- z-coordinate, [m]
        U : np.array -- velocity, [m/s]
        Re : np.array -- Reynolds number, [-]
        alpha : np.array -- angle of attack, [rad]
        alpha_c : np.array -- speed ratio, [-]
        b_c : np.array -- correlation coefficient, [-]
        p_ref : float -- reference pressure, [Pa]
        c_0 : float -- speed of sound, [m/s]
        rho_0 : float -- air density, [kg/m^3]
        cutoff : float -- radial cutoff, [-]
        probe_top : np.array -- probe location at the top surface, [-]
        probe_bot : np.array -- probe location at the bottom surface, [-]
        it : int -- number of iterations for XFoil, [-]
    """

    # Determine the Mach number
    M = U / c_0

    # Determine the Prandtl-Glauert factor
    beta = np.sqrt(1 - np.square(M))

    # Determine distance corrected for convection effects
    S_0 = np.sqrt(np.square(x) + np.square(beta) * (np.square(y) + np.square(z)))

    # Determine the angular frequency
    omega = 2 * np.pi * f

    # Determine the acoustic wavenumber
    k = omega / c_0

    # Determine the convection velocity
    U_c = U / alpha_c

    # Determine the convective wave number
    K = omega / U

    # # Determine the aerodynamic wavenumbers
    K_1 = alpha_c * K
    K_2 = k * y / S_0

    # Determine the non-dimensional wavenumbers
    K_line = K * (c/2)
    K_1_line = K_1 * (c/2)
    K_2_line = K_2 * (c/2)

    # Determine the frequency parameters
    mu_line = K_line * M / np.square(beta)
    kappa_line = np.sqrt(np.square(mu_line) - np.square(K_2_line) / np.square(beta))

    if np.any(K_2_line > K_line * M / beta):
        print("Sub-critical gust detected!")
        print("I should probably implement them :(")
        sys.exit(1)

    # Determine the correction factor
    epsilon = np.pow(1 + 1 / (4 * mu_line), -0.5)

    # Determine the average boundary layer displacement thicknesses
    delta_star_top, delta_star_bot = displacement_thickness(blade, pos, Re, M, alpha, cutoff, \
                                                            probe_top, probe_bot, it)
    delta_star = (delta_star_top + delta_star_bot) / 2

    # Determine the first-order scattering term (neglect the term np.exp(-2 * 1j * C))
    B = K_1_line + M * mu_line + kappa_line
    C = K_1_line - mu_line * (x / S_0 - M)

    f_1 = - np.exp(2 * 1j * C) / (1j * C) * ((1 + 1j) * np.exp(-2 * 1j * C) \
        * np.sqrt(B / (B - C)) * E(2 * (B - C)) - (1 + 1j) * E(2 * B) + 1)

    # Determine the second-order back-scattering correction
    f_2 = np.exp(4 * 1j * kappa_line) * (1 - (1 + 1j) * E(4 * kappa_line))
    f_2.imag = f_2.imag * epsilon

    D = kappa_line - mu_line * x / S_0

    G = (1 + epsilon) * np.exp(1j * (2 * kappa_line + D)) * np.sin(D - 2 * kappa_line) \
      / (D - 2 * kappa_line) + (1 - epsilon) * np.exp(1j * (-2 * kappa_line + D)) \
      * np.sin(D + 2 * kappa_line) / (D + 2 * kappa_line) + (1 + epsilon) * (1 - 1j) \
      / (2 * (D - 2 * kappa_line)) * np.exp(4 * 1j * kappa_line) * E(4 * kappa_line) \
      - (1 - epsilon) * (1 + 1j) / (2 * (D + 2 * kappa_line)) * np.exp(-4 * 1j * kappa_line) \
      * E(4 * kappa_line) + np.exp(2 * 1j * D) / 2 * np.sqrt(2 * kappa_line / D) * E(2 * D) \
      * ((1 + 1j) * (1 - epsilon) / (D + 2 * kappa_line) - (1 - 1j) * (1 + epsilon) \
      / (D - 2 * kappa_line))

    f_2 = f_2 - np.exp(2 * 1j * D) + 1j * (D + K_line + M * mu_line - kappa_line) * G

    # Determine the radiation integral
    Theta = np.sqrt((K_1_line + mu_line * M + kappa_line) / (K_line + mu_line * M + kappa_line))

    H = (1 + 1j) * np.exp(-4 * 1j * kappa_line) * (1 - np.square(Theta)) \
      / (2 * np.sqrt(np.pi) * (alpha_c - 1) * K_line * np.sqrt(B))

    I = f_1 + H * f_2

    # Determine the wall pressure spectrum using Schlinker's model
    Phi_pp = wall_pressure(K_1, U, delta_star, rho_0)

    # Determine the spanwise correlation length using Corcos' model
    l_y = correlation_length(omega, K_1, K_2, U_c, delta_star, b_c)

    # Determine the streamwise-integrated wavenumber spectral density of wall-pressure ﬂuctuations
    Pi_0 = (1/np.pi) * Phi_pp * l_y

    # Determine the double-sided farfield acoustic PSD. Only valid if  b >> c, AKA if we look only
    # at the priveleged oblique gust. This is investigated in (Roger and Moreau, 2009).
    S_pp = np.square(omega * z * c / (4 * np.pi * c_0 * np.square(S_0))) * 2 * np.pi * b \
         * np.square(np.abs(I)) * Pi_0

    # Determine the 1/3 octave band SPL
    SPL = 10 * np.log10(4 * np.pi * 0.232 * f * S_pp / np.square(p_ref))

    return SPL
