"""
Author:   T. Moreira da Fonte Fonseca Teles
Email:    tmoreiradafont@tudelft.nl
Date:     2025-03-14
License:  GNU GPL 3.0

Calculate the inflow noise spectra.

Classes:
    None

Functions:
    intensity
    length_scale
    directivity
    amiet
    moriarty
    inflow_noise

Exceptions:
    None
"""

import sys

import numpy as np

from misc import sears


def intensity(z, z_0, model="ZHS"):
    """
    Determine the turbulence intensity.

    Parameters:
        z : np.array -- height above the ground, [m]
        z_0 :  np.array -- surface roughness length, [m]

    Returns:
        I :  np.array -- turbulence intensity, [-]
    """

    if model == "ZHS":
        gamma = 0.24 + 0.096 * np.log10(z_0) + 0.016 * np.square(np.log10(z_0))
        I = gamma * np.log(30/z_0) / np.log(z/z_0)

    else:
        print("Unknown turbulence intensity model!")
        sys.exit(1)

    return I


def length_scale(z, z_0, model="ZHS"):
    """
    Determine the turbulence length scale.

    Parameters:
        z : np.array -- height above the ground, [m]
        z_0 :  np.array -- surface roughness length, [m]

    Returns:
        L :  np.array -- turbulence length scale, [m]
    """

    if model == "ZHS":
        L = 25 * np.power(z, 0.35) * np.pow(z_0, -0.063)

    else:
        print("Unknown turbulence length scale model!")
        sys.exit(1)

    return L


def directivity(M, theta_e, phi_e):
    """
    Determine the directivity pattern.

    Arguments:
        M : np.array -- Mach number, [-]
        theta_e : np.array -- retarded angle, [rad]
        phi_e : np.array -- retarded angle, [rad]

    Returns:
        D_L_line : np.array -- directivity array, [-]
    """

    # Determine the low-frequency directivity function
    D_L_line = np.square(np.sin(theta_e)) * np.square(np.sin(phi_e)) \
               / np.pow(1 + M * np.cos(theta_e), 4)

    return D_L_line


def amiet(f, b, c, r_e, U, alpha, L, I, c_0, rho_0):
    """
    Determine the SPL based on a flat plate.

    Arguments:
        f : np.array -- frequency, [Hz]
        b : np.array -- span, [m]
        c : np.array -- chord, [m]
        r_e : np.array -- distance, [m]
        U : np.array -- velocity, [m/s]
        alpha : np.array -- angle of attack, [rad]
        L : float -- turbulence length scale, [m]
        I : float -- turbulence intensity, [-]
        c_0 : float -- speed of sound, [m/s]
        rho_0 : float -- air density, [kg/m^3]

    Returns:
        spl_amiet : np.array -- flat plate SPL, [dB]
    """

   # Determine the high-frequency SPL
   # Use 78.4 instead of 58.4 due to units of rho_0 and c_0
    M = U / c_0
    K_x = (2*np.pi*f) / U
    K_e = 0.75 / L
    K_x_hat = K_x / K_e

    spl_h = 10 * np.log10(np.pow(M, 5) * (L*b)/(2*np.square(r_e)) \
                          * np.square(I) * np.square(rho_0) * np.pow(c_0, 4) \
                          * np.pow(K_x_hat, 3) / np.pow(1 + np.square(K_x_hat), 7/3)) + 78.4

    # Determine the low-frequency correction
    K_x_line = K_x * c / 2
    beta = np.sqrt(1 - np.square(M))

    S = sears(K_x_line / np.square(beta))

    LFC = 10 * (np.square(S) * M * np.square(K_x_line) / np.square(beta)) \
        * (1 + 9 * np.square(alpha))

    # Determine the corrected SPL
    spl_amiet = spl_h + 10 * np.log10(LFC / (1 + LFC))

    return spl_amiet


def moriarty(f, c, tc_01, tc_10, U):
    """
    Calculate the SPL correction due to the airfoil thickness.
    
    Args:
        f : np.array -- frequency, [Hz]
        c : np.array -- chord, [m]
        t_01 : np.array -- thickness at 1% c, [-]
        t_10 : np.array -- thickness at 10% c, [-]
        U : np.array -- velocity, [m/s]

    Returns:
        delta_spl : np.array -- SPL correction, [dB]
    """

    # Determine the number of frequencies and panels
    n_frequencies = f.size
    n_panels = U.size

    # Initialise the SPL matrix
    delta_spl = np.zeros((n_frequencies, n_panels))

    # Check for high Strouhal numbers
    St = (2*np.pi*f) * c / U

    if np.any(St > 75):
        print("High Strouhal numbers detected (St>75)!")

    # Determine the graph slope
    IT = tc_01 + tc_10
    SL = 1.123 * IT + 5.317 * np.square(IT)

    # Determine the SPL correction
    delta_spl = -SL * (St + 5)

    return delta_spl


def inflow_noise(f, b, c, tc_01, tc_10, x, y, z, U, alpha, c_0, rho_0, z_0, spl_correction, I_model, L_model):
    """
    Determine the inflow noise.

    Arguments:
        f : np.array -- frequency, [Hz]
        b : np.array -- span, [m]
        c : np.array -- chord, [m]
        tc_01 : np.array -- thickness at 1% c, [-]
        tc_10 : np.array -- thickness at 10% c, [-]
        x : np.array -- x-coordinate, [m]
        y : np.array -- y-coordinate, [m]
        z : np.array -- z-coordinate, [m]
        U : np.array -- velocity, [m/s]
        alpha : np.array -- angle of attack, [rad]
        c_0 : float -- speed of sound, [m/s]
        rho_0 : float -- air density, [kg/m^3]
        z_0 : float -- surface roughness length, [m]
        spl_correction : bool -- apply SPL correction?
        I_model : str -- turbulence intensity model
        L_model : str -- turbulence length scale model

    Returns:
        spl : np.array -- inflow noise SPL, [dB]
    """

    # Determine the turbulence properties
    I = intensity(z, z_0, I_model)
    L = length_scale(z, z_0, L_model)

    # Determine the retarded angles and distance
    r = np.sqrt(np.square(x) + np.square(y) + np.square(z))
    theta = np.arccos(x / r)
    phi = np.arctan2(z, y)

    theta_e = np.arccos(np.sqrt(1 - np.square(U/c_0) * np.square(np.sin(theta))) * np.cos(theta) \
                        - (U/c_0) * np.square(np.sin(theta)))
    phi_e = phi
    r_e = r * np.sin(theta) / np.sin(theta_e)

    # Determine the SPL based on a flat plate.
    spl_amiet = amiet(f, b, c, r_e, U, alpha, L, I, c_0, rho_0)

    # Determine the SPL correction
    delta_spl = moriarty(f, c, tc_01, tc_10, U)

    # Determine the directivity pattern
    D_L_line = directivity(U/c_0, theta_e, phi_e)

    # Determine the inflow noise SPL
    spl = spl_amiet + delta_spl + 10 * np.log10(D_L_line)

    # Add SPL correction
    if spl_correction:
        spl += 10

    return spl
