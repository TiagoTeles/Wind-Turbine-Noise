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


def density(p, T):
    """
    Determine the density of air.

    Arguments:
        p : float -- Pressure, [Pa]
        T : float -- Temperature, [K]

    Returns:
        rho : float -- Density, [kg/m^3]
    """

    # Determine the density
    rho = p / (287.1 * T)

    return rho


def kinematic_viscosity(p, T):
    """
    Determine the kinematic viscosity of air.

    Arguments:
        p : float -- Pressure, [Pa]
        T : float -- Temperature, [K]

    Returns:
        nu : float -- Kinematic viscosity, [m^2/s]
    """

    # Determine the dynamic viscosity
    mu = 1.458E-6 * np.power(T, 1.5) / (T + 110.4)

    # Determine the density
    rho = density(p, T)

    # Determine the kinematic viscosity
    nu = mu / rho

    return nu


def speed_of_sound(T):
    """
    Determine the speed of sound in air.

    Arguments:
        T : float -- Temperature, [K]

    Returns:
        c : float -- Speed of sound, [m/s]
    """

    # Determine the speed of sound
    c = np.sqrt(1.4 * 287.1 * T)

    return c

if __name__ == "__main__":

    # Show the density of air
    print("Rho: " + str(density(101325.0, 288.15)) + " [kg/m^3]")

    # Show the kinematic viscosity of air
    print("Nu: " + str(kinematic_viscosity(101325.0, 288.15)) + " [m^2/s]")

    # Show the speed of sound in air
    print("C: " + str(speed_of_sound(288.15)) + " [m/s]")
