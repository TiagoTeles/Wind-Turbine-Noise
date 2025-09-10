"""
Author:   T. Moreira da Fonte Fonseca Teles
Email:    tmoreiradafont@tudelft.nl
Date:     2025-09-05
License:  GNU GPL 3.0

Determine the attenuation coefficient due to atmospheric absorption.

Classes:
    None

Functions:
    attenuation_coefficient
    molar_concentration

Exceptions:
    None
"""

import matplotlib.pyplot as plt
import numpy as np

from source.atmosphere.international_standard_atmosphere import speed_of_sound
from source.constants import P_REF, T_REF, T_01, X_N, X_O, THETA_N, THETA_O
from source.settings import P_0, T_0, H_R, F_MIN, F_MAX


def attenuation_coefficient(f, p, T, h):
    """
    Determine the attenuation coefficient due to atmospheric absorption.

    Parameters:
        f : np.array -- frequency, [Hz]
        p : np.array -- pressure, [Pa]
        T : np.array -- temperature, [K]
        h : np.array -- molar concentration of water vapour, [-]

    Returns:
        alpha : np.array -- attenuation coefficient, [dB/m]
    """

    # Check for high pressures
    if np.any(p > 200000.0):
        print("High pressure detected! p > 200000.0 [Pa].")

    # Check for low temperatures
    if np.any(T < 253.15):
        print("Low temperature detected! T < 253.15 [K].")

    # Check for high temperatures
    if np.any(T > 323.15):
        print("High temperature detected! T > 323.15 [K].")

    # Check for low frequency-pressure ratios
    if np.any(f/p < 4.0E-4):
        print("Low frequency-pressure ratio detected! f/p < 4.0E-4 [Hz/Pa].")

    # Check for high frequency-pressure ratios
    if np.any(f/p > 10.0):
        print("High frequency-pressure ratio detected! f/p > 10.0 [Hz/Pa].")

    # Determine the speed of sound in air
    c_0 = speed_of_sound(T)

    # Convert the molar concentration of water vapour from [-] to [%]
    h *= 100.0

    # Determine the relaxation frequency of oxygen molecules
    f_rO = (p / P_REF) * (24.0 + 4.04E4 * h * (0.02 + h) / (0.391 + h))

    # Determine the relaxation frequency of nitrogen molecules
    f_rN = (p / P_REF) * np.pow(T / T_REF, -1.0/2.0) \
         * (9.0 + 280.0 * h * np.exp(-4.170 * (np.pow(T / T_REF, -1.0/3.0) - 1.0)))

    # Determine the maximum attenuations due to vibrational relaxation of oxygen molecules
    alpha_lambda_O = 1.559 * X_O * np.square(THETA_O / T) * np.exp(-THETA_O / T)

    # Determine the maximum attenuations due to vibrational relaxation of nitrogen molecules
    alpha_lambda_N = 1.559 * X_N * np.square(THETA_N / T) * np.exp(-THETA_N / T)

    # Determine the classical attenuation coefficient
    alpha_classical = 1.60E-10 * np.square(f) * np.pow(p / P_REF, -1.0) * np.pow(T / T_REF, 1.0/2.0)

    # Determine the attenuation coefficient due to vibrational relaxation of oxygen molecules
    alpha_vib_O = alpha_lambda_O * (f / c_0) * (2.0 * (f / f_rO) * np.pow(1.0 + np.square(f / f_rO), -1.0))

    # Determine the attenuation coefficient due to vibrational relaxation of nitrogen molecules
    alpha_vib_N = alpha_lambda_N * (f / c_0) * (2.0 * (f / f_rN) * np.pow(1.0 + np.square(f / f_rN), -1.0))

    # Determine the total attenuation coefficient
    alpha = alpha_classical + alpha_vib_O + alpha_vib_N

    return alpha


def molar_concentration(p, T, h_r):
    """
    Determine the molar concentration of water vapour.

    Parameters:
        p : np.array -- pressure, [Pa]
        T : np.array -- temperature, [K]
        h_r : np.array -- relative humidity, [-]

    Returns:
        h : np.array -- molar concentration of water vapour, [-]
    """

    # Determine the saturation pressure of water vapour
    C = -6.8346 * np.pow(T_01 / T, 1.261) + 4.6151
    p_sat = P_REF * np.pow(10.0, C)

    # Determine the molar concentration of water vapour
    h = h_r * p_sat / p

    # Check for low molar concentrations of water vapour
    if np.any(h < 5.0E-4):
        print("Low molar concentration of water vapour detected! h < 5.0E-4 [-].")

    # Check for high molar concentrations of water vapour
    if np.any(h > 0.05):
        print("High molar concentration of water vapour detected! h > 0.05 [-].")

    return h

if __name__ == "__main__":

    # Determine the molar concentration of water vapour
    h = molar_concentration(P_0, T_0, H_R)

    # Show the attenuation coefficient
    f = np.linspace(F_MIN, F_MAX, 1000)
    alpha = attenuation_coefficient(f, P_0, T_0, h)

    plt.loglog(f, alpha)
    plt.xlabel("Frequency, [Hz]")
    plt.ylabel("Attenuation Coefficient, [dB/m]")
    plt.xlim(F_MIN, F_MAX)
    plt.ylim(1.0E-5, 1.0)
    plt.grid(which="both")
    plt.show()
