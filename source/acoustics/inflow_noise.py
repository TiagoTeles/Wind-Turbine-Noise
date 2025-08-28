"""
Author:   T. Moreira da Fonte Fonseca Teles
Email:    tmoreiradafont@tudelft.nl
Date:     2025-08-28
License:  GNU GPL 3.0

Determine the inflow turbulence noise spectra.

Classes:
    None

Functions:
    flat_plate_spl
    airfoil_shape_correction
    retarded_coordinates
    inflow_noise

Exceptions:
    None
"""

import matplotlib.pyplot as plt
import numpy as np


def flat_plate_spl(f, b, c, r_e, theta_e, phi_e, U, alpha, I, L, c_0, rho_0):
    """
    Determine the flat plate SPL.

    Arguments:
        f : np.array -- frequency, [Hz]
        b : np.array -- span, [m]
        c : np.array -- chord, [m]
        r_e : np.array -- retarded distance, [m]
        theta_e : np.array -- chordwise retarded angle,
        phi_e : np.array -- spanwise retarded angle, [rad]
        U : np.array -- velocity, [m/s]
        alpha : np.array -- angle of attack, [rad]
        I : np.array -- turbulence intensity, [-]
        L : np.array -- turbulence length scale, [m]
        c_0 : float -- speed of sound, [m/s]
        rho_0 : float -- air density, [kg/m^3]

    Returns:
        spl_flat_plate : np.array -- flat plate SPL, [dB]
    """

    # Determine the Mach number
    M = U / c_0

    # Determine the angular frequency
    omega = 2 * np.pi * f

    # Determine the wavenumber in the chordwise direction
    K_x = omega / U

    # Determine the wavenumber of energy-containing eddies
    K_e = 0.75 / L

    # Determine the normalised wavenumber
    K_x_hat = K_x / K_e

    # Determine the low-frequency directivity pattern
    D_line = np.square(np.sin(theta_e)) * np.square(np.sin(phi_e)) \
               / np.pow(1 + M * np.cos(theta_e), 4)

    # Determine the high-frequency SPL (P_REF = 2E-5 [Pa])
    spl_h = 10 * np.log10(np.pow(M, 5) * (L*b)/(2*np.square(r_e)) * np.square(I) \
                          * np.square(rho_0) * np.pow(c_0, 4) * np.pow(K_x_hat, 3) \
                          / np.pow(1 + np.square(K_x_hat), 7/3) * D_line) + 78.4

    # Determine the non-dimensional wavenumber
    K_x_line = K_x * c / 2

    # Determine the Prandtl-Glauert factor
    beta = np.sqrt(1 - np.square(M))

    # Determine the compressible Sears function
    S = np.sqrt(1 / (2 * np.pi * (K_x_line / np.square(beta)) \
                     + 1 / (1 + 2.4 * (K_x_line / np.square(beta)))))

    # Determine the low-frequency correction
    LFC = 10 * np.square(S) * M * np.square(K_x_line) / np.square(beta)

    # Apply the angle of attack correction
    LFC *= 1 + 9 * np.square(alpha)

    # Determine the flat plate SPL
    spl_flat_plate = spl_h + 10 * np.log10(LFC / (1 + LFC))

    return spl_flat_plate


def airfoil_shape_correction(f, c, tc_01, tc_10, U):
    """
    Determine the airfoil shape SPL correction.
    
    Args:
        f : np.array -- frequency, [Hz]
        c : np.array -- chord, [m]
        t_01 : np.array -- thickness at x/c = 1%, [-]
        t_10 : np.array -- thickness at x/c = 10%, [-]
        U : np.array -- velocity, [m/s]

    Returns:
        delta_spl : np.array -- SPL correction, [dB]
    """

    # Determine the angular frequency
    omega = 2 * np.pi * f

    # Determine the Strouhal number
    St = omega * c / U

    # Check for high Strouhal numbers
    if np.any(St > 75):
        print("High Strouhal number detected! St > 75 [-].")

    # Determine the noise indicator
    IT = tc_01 + tc_10

    # Determine the slope parameter
    SL = 1.123 * IT + 5.317 * np.square(IT)

    # Determine the SPL correction
    delta_spl = -SL * (St + 5)

    return delta_spl


def retarded_coordinates(x, y, z, M):
    """
    Determine the retarded distance and angles.

    Arguments:
        x : np.array -- x-coordinate, [m]
        y : np.array -- y-coordinate, [m]
        z : np.array -- z-coordinate, [m]
        M : np.array -- Mach number, [-]

    Returns:
        r_e : np.array -- retarded distance, [m]
        theta_e : np.array -- chordwise retarded angle, [rad]
        phi_e : np.array -- spanwise retarded angle, [rad]
    """

    # Determine the distance
    r = np.sqrt(np.square(x) + np.square(y) + np.square(z))

    # Determine the chordwise angle
    theta = np.arccos(x / r)

    # Determine the spanwise angle
    phi = np.arctan2(z, y)

    # Determine the chordwise retarded angle
    theta_e = np.arccos(np.sqrt(1 - np.square(M) * np.square(np.sin(theta))) * np.cos(theta) \
                        - M * np.square(np.sin(theta)))

    # Determine the retarded distance
    r_e = r * np.sin(theta) / np.sin(theta_e)

    # Determine the spanwise retarded angle
    phi_e = phi

    return r_e, theta_e, phi_e


def inflow_noise(f, b, c, tc_01, tc_10, x, y, z, U, alpha, I, L, c_0, rho_0):
    """
    Determine the inflow noise SPL.

    Arguments:
        f : np.array -- frequency, [Hz]
        b : np.array -- span, [m]
        c : np.array -- chord, [m]
        tc_01 : np.array -- thickness at x/c = 1%, [-]
        tc_10 : np.array -- thickness at x/c = 10%, [-]
        x : np.array -- x-coordinate, [m]
        y : np.array -- y-coordinate, [m]
        z : np.array -- z-coordinate, [m]
        U : np.array -- velocity, [m/s]
        alpha : np.array -- angle of attack, [rad]
        I : np.array -- turbulence intensity, [-]
        L : np.array -- turbulence length scale, [m]
        c_0 : float -- speed of sound, [m/s]
        rho_0 : float -- air density, [kg/m^3]

    Returns:
        spl : np.array -- inflow noise SPL, [dB]
    """

    # Determine the Mach number
    M = U / c_0

    # Determine the retarded coordinates
    r_e, theta_e, phi_e = retarded_coordinates(x, y, z, M)

    # Determine the flat plate SPL
    spl_amiet = flat_plate_spl(f, b, c, r_e, theta_e, phi_e, U, alpha, I, L, c_0, rho_0)

    # Determine the airfoil shape SPL correction
    delta_spl = airfoil_shape_correction(f, c, tc_01, tc_10, U)

    # Determine the inflow noise SPL
    spl = spl_amiet + delta_spl

    return spl

if __name__ == "__main__":

    # Show the directivity pattern
    theta = np.linspace(0, 2*np.pi, 360)
    D_line = np.square(np.sin(theta))

    plt.polar(theta, D_line)
    plt.xlim(0, 2*np.pi)
    plt.ylim(0, 1.2)
    plt.show()
