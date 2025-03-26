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
    TODO

Exceptions:
    None
"""

import os

import numpy as np

from misc import E
from settings import XFOIL_PATH
from xfoil import XFoil


def delta_star(blade, pos, Re, M, alpha, cutoff, probe_top, probe_bot, it):
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
            top, bot = xfoil.run(path, Re[i], M[i], alpha[i], it=it)

            # Determine the displacement thickness at the probe locations
            delta_star_top[i] = np.interp(probe_top, top["x/c"], top["delta_star"])
            delta_star_bot[i] = np.interp(probe_bot, bot["x/c"], bot["delta_star"])

        else:

            # Set the displacement thickness to zero
            delta_star_top[i] = 0
            delta_star_bot[i] = 0

    return delta_star_top, delta_star_bot


def TODO(blade, pos, Re, alpha, cutoff, probe_top, probe_bot, it, f, p_ref, x, y, z, b_c,  c_0, rho_0, alpha_c, U, b, c):
    """
    TODO
    """

    # Determine the average boundary layer displacement thicknesses
    # TODO: Paper uses flat plate at aoa=0, so bl is the same on top and bot
    delta_star_top, delta_star_bot = delta_star(blade, pos, Re, M, alpha, cutoff, probe_top, probe_bot, it)
    d_star = (delta_star_top + delta_star_bot) / 2

    omega = 2 * np.pi * f
    K_x = omega / U
    M = U / c_0
    beta = np.sqrt(1 - np.square(M))
    S_0 = np.sqrt(np.square(x) + np.square(beta) * (np.square(y) + np.square(z)))
    K_2 = k * y / S_0
    U_c = U / alpha_c
    k = omega / c_0
    K = omega / U
    K_line = K * (c/2)
    mu_line = K_line * M / np.square(beta)
    K_1_line = omega  * (c/2) / U_c    # TODO: CHECK IF THIS IS CORRECT
    K_2_line = K_2 * (c/2)
    kappa_line = np.sqrt(np.square(mu_line) - np.square(K_2_line) / np.square(beta))

    if np.any(np.square(mu_line) - np.square(K_2_line) / np.square(beta) < 0):
        print("Sub-critical gust detected!")
        print("I probably need to implement them :(")

    # Deteremine the first-order scattering term
    # Neglect the term as inidcated in the paper
    epsilon = np.pow(1 + 1 / (4 * mu_line), -0.5)

    B = K_1_line + M * mu_line + kappa_line
    C = K_1_line - mu_line * (x / S_0 - M)

    uncorrected_1 = (1 + 1j) * np.exp(-2 * 1j * C) * np.sqrt(B/(B-C)) * E(2 * (B - C)) - (1 + 1j) * E(2 * B) + 1
    corrected_1 = np.zeros(uncorrected_1.shape, dtype=complex)
    corrected_1.real = uncorrected_1.real
    corrected_1.imag = uncorrected_1.imag * epsilon

    f_1 = - np.exp(2 * 1j * C) / (1j * C) * corrected_1

    # Determine the second-order back-scattering correction
    uncorrected_2 = np.exp(4 * 1j * kappa_line) * (1 - (1 + 1j) * E(4 * kappa_line))

    corrected_2 = np.zeros(uncorrected_2.shape, dtype=complex)
    corrected_2.real = uncorrected_2.real
    corrected_2.imag = uncorrected_2.imag * epsilon

    D = kappa_line - mu_line * x / S_0

    G = (1 + epsilon) * np.exp(1j * (2 * kappa_line + D)) * np.sin(D - 2 * kappa_line) / (D - 2 * kappa_line) + (1 - epsilon) * np.exp(1j * (-2 * kappa_line + D)) * np.sin(D + 2*kappa_line) / (D + 2 * kappa_line)
    G += (1 + epsilon) * (1 - 1j) / (2 * (D - 2 * kappa_line)) * np.exp(4 * 1j * kappa_line) * E(4 * kappa_line) - (1 - epsilon) * (1 + 1j) / (2 * (D + 2 * kappa_line)) * np.exp(-4 * 1j * kappa_line) * E(4 * kappa_line)
    G += np.exp(2 * 1j * D) / 2 * np.sqrt(2 * kappa_line / D) * E(2 * D) * ((1 + 1j) * (1 - epsilon) / (D + 2 * kappa_line) - (1 - 1j) * (1 + epsilon) / (D - 2 * kappa_line))

    f_2 = corrected_2 - np.exp(2 * 1j * D) + 1j * (D + K_line + M * mu_line - kappa_line) * G

    # Determine the radiation integral
    Theta = np.sqrt((K_1_line + M * mu_line + kappa_line) / (K_line + M * mu_line + kappa_line))
    H = (1 + 1j) * np.exp(-4 * 1j * kappa_line) * (1 - np.square(Theta)) / (2 * np.sqrt(np.pi) * (alpha_c - 1) * K_line * np.sqrt(B))
    I = f_1 + H * f_2

    # Determine the wall pressure spectrum using an empirical model from data from an airfoil
    # TODO: Is there a newer model for F or Phi_pp
    omega_tilde = K_x * d_star
    F = (33.28 * omega_tilde) / (1 - 5.489 * omega_tilde + 36.74 * np.square(omega_tilde) + 0.1505 * np.pow(omega_tilde, 5))
    Phi_pp = np.square(0.5 * rho_0 * np.square(U)) * (d_star / U) * 2E-5 * F

    # Determine the spanwise correlation length using Corcos' model
    # TODO: Is this model good enough?
    l_y = (omega / (b_c * U_c)) / (np.square(K_2) + np.square(omega) / np.square(b_c * U_c))
    l_y *= (1 - np.exp(-np.square(omega_tilde) / 0.09)) # Prof. Dr. D. Casalino correction term

    # Determine the streamwise-integrated wavenumber spectral density of wall-pressure ﬂuctuations
    Pi_0 = (1/np.pi) * Phi_pp * l_y

    # Determine the double-sided farfield acoustic PSD. This equation simplifies if  b >> c, AKA if
    # we look only at the priveleged oblique gust. This is investigated in the validation of Part 2
    # TODO: Does the factor of 2 before np.pi model both airfoil surfaces? Comparing equations seems
    # to indicate that it does. Also stated in Roger and Moreau (2005)
    S_pp = np.square(omega * z * c / (4 * np.pi * c_0 * np.square(S_0))) * 2 * np.pi * b \
         * np.square(np.abs(I)) * Pi_0

    # Determine the 1/3 octave band SPL
    SPL = 10 * np.log10((2*np.pi*(0.232*f)) * (2*S_pp) / np.square(p_ref))

    return SPL
