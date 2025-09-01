"""
Author:   T. Moreira da Fonte Fonseca Teles
Email:    tmoreiradafont@tudelft.nl
Date:     2025-08-25
License:  GNU GPL 3.0

Functions from the International Standard Atmosphere.

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

from source.constants import R, GAMMA, S, BETA_S
from source.settings import P_0, T_0


def density(p, T):
    """
    Determine the density of air.

    Arguments:
        p : float -- pressure, [Pa]
        T : float -- temperature, [K]

    Returns:
        rho : float -- density, [kg/m^3]
    """

    # Determine the density
    rho = p / (R * T)

    return rho


def kinematic_viscosity(p, T):
    """
    Determine the kinematic viscosity of air.

    Arguments:
        p : float -- pressure, [Pa]
        T : float -- temperature, [K]

    Returns:
        nu : float -- kinematic viscosity, [m^2/s]
    """

    # Determine the dynamic viscosity
    mu = BETA_S * np.power(T, 3/2) / (T + S)

    # Determine the density
    rho = density(p, T)

    # Determine the kinematic viscosity
    nu = mu / rho

    return nu


def speed_of_sound(T):
    """
    Determine the speed of sound in air.

    Arguments:
        T : float -- temperature, [K]

    Returns:
        c : float -- speed of sound, [m/s]
    """

    # Determine the speed of sound
    c = np.sqrt(GAMMA * R * T)

    return c

if __name__ == "__main__":

    # Show the density of air
    print("Rho: " + str(density(P_0, T_0)) + " [kg/m^3]")

    # Show the kinematic viscosity of air
    print("Nu: " + str(kinematic_viscosity(P_0, T_0)) + " [m^2/s]")

    # Show the speed of sound in air
    print("C: " + str(speed_of_sound(T_0)) + " [m/s]")
