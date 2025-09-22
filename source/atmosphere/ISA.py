"""
Author:   T. Moreira da Fonte Fonseca Teles
Email:    tmoreiradafont@tudelft.nl
Date:     2025-09-17
License:  GNU GPL 3.0

Determine the international standard atmosphere.

Classes:
    None

Functions:
    density
    kinematic_viscosity
    speed_of_sound

Exceptions:
    None
"""

import numpy as np

from source.constants import R, S, BETA_S, GAMMA
from source.settings import P_0, T_0


def density(p, T):
    """
    Determine the density of air.

    Parameters:
        p : np.ndarray -- pressure, [Pa]
        T : np.ndarray -- temperature, [K]

    Returns:
        rho : np.ndarray -- density, [kg/m^3]
    """

    # Determine the density
    rho = p / (R * T)

    return rho


def kinematic_viscosity(T, rho):
    """
    Determine the kinematic viscosity of air.

    Parameters:
        T : np.ndarray -- temperature, [K]
        rho : np.ndarray -- density, [kg/m^3]

    Returns:
        nu : np.ndarray -- kinematic viscosity, [m^2/s]
    """

    # Determine the dynamic viscosity
    mu = BETA_S * np.power(T, 3.0 / 2.0) / (T + S)

    # Determine the kinematic viscosity
    nu = mu / rho

    return nu


def speed_of_sound(T):
    """
    Determine the speed of sound in air.

    Parameters:
        T : np.ndarray -- temperature, [K]

    Returns:
        c : np.ndarray -- speed of sound, [m/s]
    """

    # Determine the speed of sound
    c = np.sqrt(GAMMA * R * T)

    return c

if __name__ == "__main__":

    # Show the density of air
    rho = density(P_0, T_0)
    print("Rho: " + str(rho) + " [kg/m^3]")

    # Show the kinematic viscosity of air
    nu = kinematic_viscosity(T_0, rho)
    print("Nu: " + str(nu) + " [m^2/s]")

    # Show the speed of sound in air
    c = speed_of_sound(T_0)
    print("C: " + str(c) + " [m/s]")
