"""
Author:   T. Moreira da Fonte Fonseca Teles
Email:    tmoreiradafont@tudelft.nl
Date:     2025-10-28
License:  GNU GPL 3.0

Determine the inflow turbulence noise spectra.

Classes:
    None

Functions:
    flat_plate_spl
    airfoil_shape_correction
    leading_edge_noise

Exceptions:
    None
"""

import matplotlib.pyplot as plt
import numpy as np


def flat_plate_spl(f, b, c, I, L, U, alpha, x, y, z, c_0, rho_0):
    """
    Determine the flat plate SPL.

    Parameters:
        f : np.ndarray -- frequency, [Hz]
        b : np.ndarray -- span, [m]
        c : np.ndarray -- chord, [m]
        I : np.ndarray -- turbulence intensity, [-]
        L : np.ndarray -- turbulence length scale, [m]
        U : np.ndarray -- velocity, [m/s]
        alpha : np.ndarray -- angle of attack, [rad]
        x : np.ndarray -- x coordinate, [m]
        y : np.ndarray -- y coordinate, [m]
        z : np.ndarray -- z coordinate, [m]
        c_0 : float -- speed of sound, [m/s]
        rho_0 : float -- density, [kg/m^3]

    Returns:
        spl_flat_plate : np.ndarray -- flat plate SPL, [dB]
    """

    # Determine the Mach number
    M = U / c_0

    # Determine the retarded distance
    r = np.sqrt(np.square(x) + np.square(y) + np.square(z))

    # Determine the angular frequency
    omega = 2.0 * np.pi * f

    # Determine the wavenumber in the chordwise direction
    K_x = omega / U

    # Determine the wavenumber of energy-containing eddies
    K_e = 0.75 / L

    # Determine the normalised wavenumber
    K_x_hat = K_x / K_e

    # Determine the retarded chordwise angle
    theta = np.arccos(x / r)

    # Determine the retarded spanwise angle
    phi = np.arctan2(z, y)

    # Determine the low-frequency directivity pattern
    D_line = np.square(np.sin(theta)) * np.square(np.sin(phi)) \
           / np.pow(1.0 + M * np.cos(theta), 4.0)

    # Determine the high-frequency SPL (P_REF_AIR = 2.0E-5 [Pa])
    spl_h = 10.0 * np.log10(np.pow(M, 5.0) * (L * b) / (2.0 * np.square(r)) * np.square(I) \
                            * np.square(rho_0) * np.pow(c_0, 4.0) * np.pow(K_x_hat, 3.0) \
                            / np.pow(1.0 + np.square(K_x_hat), 7.0 / 3.0) * D_line) + 78.4

    # Determine the non-dimensional wavenumber
    K_x_line = K_x * c / 2.0

    # Determine the Prandtl-Glauert factor
    beta = np.sqrt(1.0 - np.square(M))

    # Determine the compressible Sears function
    S = np.sqrt(1.0 / (2.0 * np.pi * (K_x_line / np.square(beta)) + 1.0 \
                       / (1.0 + 2.4 * (K_x_line / np.square(beta)))))

    # Determine the low-frequency correction
    LFC = 10.0 * np.square(S) * M * np.square(K_x_line) / np.square(beta)

    # Apply the angle of attack correction
    LFC *= 1.0 + 9.0 * np.square(alpha)

    # Determine the flat plate SPL
    spl_flat_plate = spl_h + 10.0 * np.log10(LFC / (1.0 + LFC))

    return spl_flat_plate


def airfoil_shape_correction(f, c, t_01, t_10, U):
    """
    Determine the airfoil shape SPL correction.

    Args:
        f : np.ndarray -- frequency, [Hz]
        c : np.ndarray -- chord, [m]
        t_01 : np.ndarray -- thickness at x/c = 0.01, [-]
        t_10 : np.ndarray -- thickness at x/c = 0.10, [-]
        U : np.ndarray -- velocity, [m/s]

    Returns:
        delta_spl : np.ndarray -- SPL correction, [dB]
    """

    # Determine the angular frequency
    omega = 2.0 * np.pi * f

    # Determine the Strouhal number
    St = omega * c / U

    # Check for high Strouhal numbers
    if np.any(St > 75.0):
        print("High Strouhal number detected! St > 75 [-].")

    # Determine the noise indicator
    IT = t_01 + t_10

    # Determine the slope parameter
    SL = 1.123 * IT + 5.317 * np.square(IT)

    # Determine the SPL correction
    delta_spl = -SL * (St + 5.0)

    return delta_spl


def leading_edge_noise(f, b, c, I, L, t_01, t_10, U, alpha, x, y, z, c_0, rho_0):
    """
    Determine the inflow turbulence noise SPL.

    Parameters:
        f : np.ndarray -- frequency, [Hz]
        b : np.ndarray -- span, [m]
        c : np.ndarray -- chord, [m]
        I : np.ndarray -- turbulence intensity, [-]
        L : np.ndarray -- turbulence length scale, [m]
        t_01 : np.ndarray -- thickness at x/c = 0.01, [-]
        t_10 : np.ndarray -- thickness at x/c = 0.10, [-]
        U : np.ndarray -- velocity, [m/s]
        alpha : np.ndarray -- angle of attack, [rad]
        x : np.ndarray -- x coordinate, [m]
        y : np.ndarray -- y coordinate, [m]
        z : np.ndarray -- z coordinate, [m]
        c_0 : float -- speed of sound, [m/s]
        rho_0 : float -- air density, [kg/m^3]

    Returns:
        spl : np.ndarray -- inflow turbulence noise SPL, [dB]
    """

    # Determine the flat plate SPL
    spl_amiet = flat_plate_spl(f, b, c, I, L, U, alpha, x, y, z, c_0, rho_0)

    # Determine the airfoil shape SPL correction
    delta_spl = airfoil_shape_correction(f, c, t_01, t_10, U)

    # Determine the inflow turbulence noise SPL
    spl = spl_amiet + delta_spl

    return spl

if __name__ == "__main__":

    # Show the directivity pattern
    theta = np.linspace(0.0, 2.0 * np.pi, 360)
    D_line = np.square(np.sin(theta)) / np.pow(1.0 + 0.3 * np.cos(theta), 4.0)

    plt.polar(theta, D_line)
    plt.xlim(0.0, 2.0 * np.pi)
    plt.ylim(0.0, 1.6)
    plt.show()
