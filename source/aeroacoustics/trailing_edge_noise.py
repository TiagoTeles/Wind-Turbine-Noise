"""
Author:   T. Moreira da Fonte Fonseca Teles
Email:    tmoreiradafont@tudelft.nl
Date:     2025-02-19
License:  GNU GPL 3.0

Calculate the TBLTE noise.

Classes:
    None

Functions:
    TBD

Exceptions:
    None
"""

import os

import numpy as np

from misc import E
from xfoil import XFoil


def trailing_edge_noise(blade, f, x, y, z, U, Re, alpha, b, c, c_0, rho_0, p_ref, xfoil_path, probe_top, probe_bot):
    
    # Determine the number of frequencies and panels
    n_panels = U.size

    delta_f = 0.232 * f

    omega = 2 * np.pi * f

    M = U / c_0

    beta = np.sqrt(1 - np.square(M))

    S_0 = np.sqrt(np.square(x) + np.square(beta) * (np.square(y) + np.square(z)))

    K_x = omega / U

    delta_star_top = np.zeros(n_panels)[np.newaxis, :, np.newaxis]
    delta_star_bot = np.zeros(n_panels)[np.newaxis, :, np.newaxis]

    # Iterate through each panel
    for i in range(n_panels):

        print(i)

        airfoil = blade.data["airfoil"].iloc[i]

        # TODO: Is suction side always the top side? Is pressure side always the bottom side?
        if "Circular" not in airfoil.attributes["AIRFOILNAME"] and airfoil.thickness(0.5) < 0.3:
            xfoil = XFoil(xfoil_path, cwd=os.path.dirname(blade.path))

            path = os.path.relpath(airfoil.path, os.path.dirname(blade.path))
            top, bot = xfoil.run(path, float(Re[0, i, 0]), float(M[0, i, 0]), float(alpha[0, i, 0]), it=1000)

            delta_star_top[:, i, :] = np.interp(probe_top, top["x/c"], top["delta_star"])
            delta_star_bot[:, i, :] = np.interp(probe_bot, bot["x/c"], bot["delta_star"])

        else:
            # delta_star_top[:, i, :] = np.nan
            # delta_star_bot[:, i, :] = np.nan
            delta_star_top[:, i, :] = 0
            delta_star_bot[:, i, :] = 0

    # TODO: For now average the delta_star_top and delta_star_bot values, and use factor of 4 in SPL calculation
    # TODO: This is because for a flat plate delta_star_top == delta_star_bot

    delta_star = delta_star_top + delta_star_bot

    omega_tilde = K_x * delta_star

    # TODO: IS THIS THE BEST WAY TO CALCULATE F? IS THERE A NEWER MODEL
    F = (33.28 * omega_tilde) / (1 - 5.489 * omega_tilde + 36.74 * np.square(omega_tilde) + 0.1505 * np.pow(omega_tilde, 5))

    # TODO: ANY NEWER EMPITICAL FORMULATION? IS 2E-5 related to P_REF???
    Phi_pp = np.square(0.5 * rho_0 * np.square(U)) * (delta_star / U) * 2E-5 * F

    # plt.plot(f[:, 0, 0], 10 * np.log10(Phi_pp[:, 15, 0]/ np.square(2E-5)))
    # plt.xscale("log")
    # plt.yscale("linear")
    # plt.grid(which="both")
    # plt.show()

    # TODO: FIND VALUE OF b_c. Roger claims it is a constant, but does not provide a value
    b_c = 1

    # TODO: FIND VALUE OF alpha_c. Roger does not provide a value
    alpha_c = 1 / 0.8

    U_c = U / alpha_c

    k = omega / c_0

    K_2 = k * y / S_0

    l_y = (omega / (b_c * U_c)) / (np.square(K_2) + np.square(omega) / np.square(b_c * U_c))

    # Prof. Casalino correction term
    l_y *= (1 - np.exp(-np.square(omega_tilde) / 0.09))

    Pi_0 = (1/np.pi) * Phi_pp * l_y

    # TODO: Why not U_c
    K = omega / U

    K_line = K * (c/2)

    mu_line = K_line * M / np.square(beta)

    epsilon = np.pow(1 + 1 / (4 * mu_line), -0.5)

    K_2_line = K_2 * (c/2)

    kappa_line = np.sqrt(np.square(mu_line) - np.square(K_2_line) / np.square(beta))

    if np.any(np.square(mu_line) - np.square(K_2_line) / np.square(beta) < 0):
        print("Sub-critical gust detected!")
        print("I probably need to implement them :(")


    K_1_line = K_line * alpha_c

    B = K_1_line + M * mu_line + kappa_line
    C = K_1_line - mu_line * (x / S_0 - M)

    # Determine supercritical gust contribution
    # Remove  "-np.exp(-2*1j*C)" as indicated in paper
    # TODO: DOES THIS PART NEED TO BE CORRECTED? HAS {} but not ^c
    uncorrected_1 = (1 + 1j) * np.exp(-2 * 1j * C) * np.sqrt(B/(B-C)) * E(2 * (B - C)) - (1 + 1j) * E(2*B) +1

    corrected_1 = np.zeros(uncorrected_1.shape, dtype=complex)
    corrected_1.real = uncorrected_1.real
    corrected_1.imag = uncorrected_1.imag * epsilon

    part_1 = - np.exp(2 * 1j * C) / (1j * C) * corrected_1

    uncorrected_2 = np.exp(4 * 1j * kappa_line) * (1 - (1 + 1j) * E(4 * kappa_line))

    corrected_2 = np.zeros(uncorrected_2.shape, dtype=complex)
    corrected_2.real = uncorrected_2.real
    corrected_2.imag = uncorrected_2.imag * epsilon

    D = kappa_line - mu_line * x / S_0

    G = (1 + epsilon) * np.exp(1j * (2 * kappa_line + D)) * np.sin(D - 2 * kappa_line) / (D - 2 * kappa_line) + (1 - epsilon) * np.exp(1j * (-2 * kappa_line + D)) * np.sin(D + 2*kappa_line) / (D + 2 * kappa_line)
    G += (1 + epsilon) * (1 - 1j) / (2 * (D - 2 * kappa_line)) * np.exp(4 * 1j * kappa_line) * E(4 * kappa_line) - (1 - epsilon) * (1 + 1j) / (2 * (D + 2 * kappa_line)) * np.exp(-4 * 1j * kappa_line) * E(4 * kappa_line)
    G += np.exp(2 * 1j * D) / 2 * np.sqrt(2 * kappa_line / D) * E(2 * D) * ((1 + 1j) * (1 - epsilon) / (D + 2 * kappa_line) - (1 - 1j) * (1 + epsilon) / (D - 2 * kappa_line))

    part_2 = corrected_2 - np.exp(2 * 1j * D) + 1j * (D + K_line + M * mu_line - kappa_line) * G

    Theta = np.sqrt((K_1_line + M * mu_line + kappa_line) / (K_line + M * mu_line + kappa_line))

    H = (1 + 1j) * np.exp(-4 * 1j * kappa_line) * (1 - np.square(Theta)) / (2 * np.sqrt(np.pi) * (alpha_c - 1) * K_line * np.sqrt(B))

    I = part_1 + H * part_2

    # TODO: Check if 2 or 4
    S_pp = np.square(omega * z * c / (4 * np.pi * c_0 * np.square(S_0))) \
           * 2 * np.pi * b * np.square(np.abs(I)) * Pi_0

    SPL = 10 * np.log10(4 * np.pi * delta_f * S_pp / np.square(p_ref))

    # return SPL
    return SPL
