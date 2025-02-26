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

import numpy as np


def integral():
    pass

def roger_moreau():

    f = 500
    x = 1
    y = 1
    z = 1
    b = 1
    C_0 = 340
    RHO_0 = 1.225
    L = 1
    U = 1
    V = 1
    c = 1
    epsilon = 1
    alpha = 1/0.8
    b_c = 1
    delta_star = 1

    p_ref = 2E-5

    omega_tilde = omega * delta_star / U

    F = (33.28 * omega_tilde)/(1 - 5.489 * omega_tilde + 36.74 * np.square(omega_tilde) + 0.1505*np.pow(omega_tilde, 5))

    Phi_pp = np.square(0.5 * RHO_0 * np.square(U)) * (delta_star / U) * 2E-5 * F

    l_y = (omega / (b_c * U_c)) / (np.square(K_2) + np.square(omega) / np.square(b_c * U_c))

    Pi_0 = (1/np.pi) * Phi_pp * l_y

    omega = 2 * np.pi * f

    M = U / C_0

    beta = np.sqrt(1 - np.square(M))

    S_0 = np.sqrt(np.square(x) + np.square(beta) * (np.square(y) + np.square(z)))

    K_1_line = alpha * K_line

    U_c = U / alpha

    K = omega / U_c 

    K_line = K * b

    mu_line = K_line * M / np.square(beta)

    C = K_1_line - mu_line * (x / S_0 - M)

    B = K_1_line + M * mu_line + kappa_line

    Theta = np.sqrt((K_1_line + M * mu_line + kappa_line)/ (K_line + M * mu_line + kappa_line))

    H = (1 + 1j) * np.exp(-4 * 1j * kappa_line) * (1 - np.square(Theta)) / (2 * np.sqrt(np.pi) * (alpha - 1) * K_line * np.sqrt(B))

    def E():
        return 0
    
    k = omega / C_0

    K_2 = k * y / S_0     # TODO: CHECK THIS

    K_2_line = K_2 * b

    kappa_line = np.sqrt(np.square(mu_line) - np.square(K_2_line) / np.square(beta))

    D = kappa_line - mu_line * x / S_0

    G = (1 + epsilon) * np.exp(1j * (2 * kappa_line + D)) * (np.sin(D - 2 * kappa_line))/(D - 2 * kappa_line) + (1 - epsilon) * np.exp(1j * (-2*kappa_line + D)) * (np.sin(D+2*kappa_line))/(D + 2 * kappa_line)
    G += (1 + epsilon) * (1 - 1j) / (2 * (D - 2 * kappa_line)) * np.exp(4 * 1j * kappa_line) * E(4 * kappa_line) - (1-epsilon) * (1+epsilon) / (2*(D + 2 * kappa_line)) * np.exp(-4 * 1j * kappa_line) * E(4 * kappa_line)
    G += np.exp(2 * 1j * D) / 2 * np.sqrt(2*kappa_line/D) * E(2*D) * (((1+1j)*(1-epsilon))/(D+2*kappa_line) - ((1-1j)*(1+1j))/(D-2*kappa_line))

    part_1 = - np.exp(2 * 1j * C) / (1j * C) * ((1 + 1j) * np.exp(-2 * 1j * C) * np.sqrt(B/(B-C)) * E(2 * (B - C)) - (1 + 1j) * E(2*B) +1 - np.exp(-2*1j*C))

    part_2 = np.imag(np.exp(4 * 1j * kappa_line) * (1 - (1+ 1j) * E(4 * kappa_line))) * epsilon - np.exp(2 * 1j * D) + 1j * (D + K_line + M * mu_line - kappa_line) * G

    
    I = part_1 + H * part_2

    S_pp = np.square(omega * z * b / (2* np.pi * C_0 * np.square(S_0))) \
           * 2 * np.pi * L * np.square(I) * Pi_0
    
    delta_f = 0.232 * f

    SPL = 10 * np.log10(4 * np.pi * delta_f * S_pp / np.square(p_ref))

    return SPL