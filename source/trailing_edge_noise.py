"""
Author:   T. Moreira da Fonte Fonseca Teles
Email:    tmoreiradafont@tudelft.nl
Date:     2025-03-18
License:  GNU GPL 3.0

Calculate the TBLTE noise.

Classes:
    None

Functions:
    TBD TODO

Exceptions:
    None
"""

import os

import numpy as np

from misc import E
from settings import XFOIL_PATH
from xfoil import XFoil


def delta_star(blade, Re, M, alpha, probe_top, probe_bot, cutoff, it):
    """
    Determine the boundary layer displacement thickness.

    Arguments:
        blade : Blade -- blade object
        Re : np.array -- Reynolds number, [-]
        M : np.array -- Mach number, [-]
        alpha : np.array -- angle of attack, [rad]
        probe_top : np.array -- probe location at the top surface, [-]
        probe_bot : np.array -- probe location at the bottom surface, [-]
        cutoff : float -- thickness cutoff for the airfoil, [-]
        it : int -- number of iterations for XFoil

    Returns:
        delta_star_top : np.array -- boundary layer displacement thickness on the top side, [m]
        delta_star_bot : np.array -- boundary layer displacement thickness on the bottom side, [m]
    """

    delta_star_top = np.zeros(Re.shape)
    delta_star_bot = np.zeros(Re.shape)

    # Iterate through each panel
    for i in range(Re.shape[1]):

        airfoil = blade.data["airfoil"].iloc[i]

        if airfoil.thickness(0.5) < cutoff:
            
            # Determine the path to the airfoil
            cwd = os.path.dirname(blade.path)
            path = os.path.relpath(airfoil.path, cwd)

            # Determine Re, M, and alpha
            Re = float(Re[0, i, 0])
            M = float(M[0, i, 0])
            alpha = float(alpha[0, i, 0])

            # Run XFoil
            xfoil = XFoil(XFOIL_PATH, cwd)
            top, bot = xfoil.run(path, Re, M, alpha, it=it)

            # Determine the displacement thickness at the probe locations
            delta_star_top[:, i, :] = np.interp(probe_top, top["x/c"], top["delta_star"])
            delta_star_bot[:, i, :] = np.interp(probe_bot, bot["x/c"], bot["delta_star"])

        else:

            # Set the displacement thickness to zero
            delta_star_top[:, i, :] = 0
            delta_star_bot[:, i, :] = 0

    return delta_star_top, delta_star_bot




def f1(K_1_line, mu_line, kappa_line, M, x, S_0, epsilon):
    """
    Determine Eq. 13 from (Roger and Moreau, 2005). Remove the "-np.exp(-2*1j*C)" 
    term as indicated in (Casalino et al., 2022).

    Arguments:
        K_1_line : np.array -- non-dimensional wavenumber, [-]
        mu_line : np.array -- non-dimensional wavenumber, [-]
        kappa_line : np.array -- non-dimensional wavenumber, [-]
        M : np.array -- Mach number, [-]
        x : np.array -- x-coordinate, [m]
        S_0 : np.array -- corrected distance for convection effects, [m]
        epsilon : float -- correction factor, [-]
    """

    B = K_1_line + M * mu_line + kappa_line
    C = K_1_line - mu_line * (x / S_0 - M)

    uncorrected_1 = (1 + 1j) * np.exp(-2 * 1j * C) * np.sqrt(B/(B-C)) * E(2 * (B - C)) - (1 + 1j) * E(2 * B) + 1

    corrected_1 = np.zeros(uncorrected_1.shape, dtype=complex)
    corrected_1.real = uncorrected_1.real
    corrected_1.imag = uncorrected_1.imag * epsilon

    return - np.exp(2 * 1j * C) / (1j * C) * corrected_1


def f2(K_line, mu_line, kappa_line, M, x, S_0, epsilon):
    """
    Determine Eq. 14 from (Roger and Moreau, 2005).

    Arguments:
        K_1_line : np.array -- non-dimensional wavenumber, [-]
        mu_line : np.array -- non-dimensional wavenumber, [-]
        kappa_line : np.array -- non-dimensional wavenumber, [-]
        M : np.array -- Mach number, [-]
        x : np.array -- x-coordinate, [m]
        S_0 : np.array -- corrected distance for convection effects, [m]
        epsilon : float -- correction factor, [-]
    """

    uncorrected_2 = np.exp(4 * 1j * kappa_line) * (1 - (1 + 1j) * E(4 * kappa_line))

    corrected_2 = np.zeros(uncorrected_2.shape, dtype=complex)
    corrected_2.real = uncorrected_2.real
    corrected_2.imag = uncorrected_2.imag * epsilon

    D = kappa_line - mu_line * x / S_0

    G = (1 + epsilon) * np.exp(1j * (2 * kappa_line + D)) * np.sin(D - 2 * kappa_line) / (D - 2 * kappa_line) + (1 - epsilon) * np.exp(1j * (-2 * kappa_line + D)) * np.sin(D + 2*kappa_line) / (D + 2 * kappa_line)
    G += (1 + epsilon) * (1 - 1j) / (2 * (D - 2 * kappa_line)) * np.exp(4 * 1j * kappa_line) * E(4 * kappa_line) - (1 - epsilon) * (1 + 1j) / (2 * (D + 2 * kappa_line)) * np.exp(-4 * 1j * kappa_line) * E(4 * kappa_line)
    G += np.exp(2 * 1j * D) / 2 * np.sqrt(2 * kappa_line / D) * E(2 * D) * ((1 + 1j) * (1 - epsilon) / (D + 2 * kappa_line) - (1 - 1j) * (1 + epsilon) / (D - 2 * kappa_line))

    return corrected_2 - np.exp(2 * 1j * D) + 1j * (D + K_line + M * mu_line - kappa_line) * G

def farfield_PSD(omega, b, c, z, S_0, c_0, I, Pi_0):
    """
    Determine the farfield PSD.

    Arguments:
        omega : np.array -- angular frequency, [rad/s]
        b : np.array -- blade span, [m]
        c : np.array -- chord length, [m]
        z : np.array -- z-coordinate, [m]
        S_0 : np.array -- corrected distance for convection effects, [m]
        c_0 : float -- speed of sound, [m/s]
        I : np.array -- radiation integral, [TODO]
        Pi_0 : np.array -- streamwise-integrated wavenumber of spectral 
                           density of wall-pressure ﬂuctuations, [TODO]

    Returns:
        S_pp : np.array -- farfield PSD, [TODO]
    """

    # TODO: Check if 2 or 4
    S_pp = np.square(omega * z * c / (2 * np.pi * c_0 * np.square(S_0))) * 2 * np.pi * b \
         * np.square(np.abs(I)) * Pi_0

    return S_pp


def trailing_edge_noise(blade, f, x, y, z, U, Re, alpha, b, c, c_0, rho_0, p_ref, probe_top, probe_bot, cutoff, it):
    """
    Determine the trailing edge noise.

    Arguments:
        blade : Blade -- blade object
        f : float -- frequency, [Hz]
        x : np.array -- x-coordinate, [m]
        y : np.array -- y-coordinate, [m]
        z : np.array -- z-coordinate, [m]
        U : float -- free-stream velocity, [m/s]
        Re : np.array -- Reynolds number, [-]
        alpha : np.array -- angle of attack, [rad]
        b : float -- blade span, [m]
        c : float -- chord length, [m]
        c_0 : float -- speed of sound, [m/s]
        rho_0 : float -- air density, [kg/m^3]
        p_ref : float -- reference pressure, [Pa]
        probe_top : float -- probe location at the top surface, [-]
        probe_bot : float -- probe location at the bottom surface, [-]
        cutoff : float -- thickness cutoff for the airfoil, [-]
        it : int -- number of iterations for XFoil

    Returns:
        SPL : np.array -- TBLTE noise SPL, [dB]
    """

    M = U / c_0

    # Determine the boundary layer displacement thicknesses
    delta_star_top, delta_star_bot = delta_star(blade, Re, M, alpha, probe_top, probe_bot, cutoff, it)
    delta_star = (delta_star_top + delta_star_bot)/2

    # Determine F (TODO: IS THIS THE BEST WAY TO CALCULATE F? IS THERE A NEWER MODEL)
    omega = 2 * np.pi * f
    K_x = omega / U
    omega_tilde = K_x * delta_star

    F = (33.28 * omega_tilde) / (1 - 5.489 * omega_tilde + 36.74 * np.square(omega_tilde) + 0.1505 * np.pow(omega_tilde, 5))

    # Determine the wall pressure spectrum
    Phi_pp = np.square(0.5 * rho_0 * np.square(U)) * (delta_star / U) * 2E-5 * F

    # Determine the radiation integral
    K = omega / U                                       # TODO: Why not U_c
    beta = np.sqrt(1 - np.square(M))
    K_line = K * (c/2)
    k = omega / c_0
    K_2 = k * y / S_0
    mu_line = K_line * M / np.square(beta)
    alpha_c = 1 / 0.8
    K_1_line = K_line * alpha_c
    K_2_line = K_2 * (c/2)
    kappa_line = np.sqrt(np.square(mu_line) - np.square(K_2_line) / np.square(beta))
    S_0 = np.sqrt(np.square(x) + np.square(beta) * (np.square(y) + np.square(z)))
    epsilon = np.pow(1 + 1 / (4 * mu_line), -0.5)

    if np.any(np.square(mu_line) - np.square(K_2_line) / np.square(beta) < 0):
        print("Sub-critical gust detected!")
        print("I probably need to implement them :(")

    f_1 = f1(K_1_line, mu_line, kappa_line, M, x, S_0, epsilon)
    f_2 = f2(K_line, mu_line, kappa_line, M, x, S_0, epsilon)

    B = K_1_line + M * mu_line + kappa_line
    Theta = np.sqrt((K_1_line + M * mu_line + kappa_line) / (K_line + M * mu_line + kappa_line))
    H = (1 + 1j) * np.exp(-4 * 1j * kappa_line) * (1 - np.square(Theta)) / (2 * np.sqrt(np.pi) * (alpha_c - 1) * K_line * np.sqrt(B))

    I = f_1 + H * f_2

    # Determine the correlation length
    b_c = 1
    U_c = U / alpha_c
    l_y = (omega / (b_c * U_c)) / (np.square(K_2) + np.square(omega) / np.square(b_c * U_c))
    l_y *= (1 - np.exp(-np.square(omega_tilde) / 0.09)) # Prof. Casalino correction term

    # Determine the streamwise-integrated wavenumber of spectral density of wall-pressure ﬂuctuations
    Pi_0 = (1/np.pi) * Phi_pp * l_y

    # Determine the farfield PSD
    S_pp = farfield_PSD(omega, b, c, z, S_0, c_0, I, Pi_0)

    # Determine the farfield SPL
    delta_f = 0.232 * f
    SPL = 10 * np.log10(4 * np.pi * delta_f * S_pp / np.square(p_ref))

    return SPL
