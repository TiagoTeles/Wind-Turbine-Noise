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


def density(p, T, R):
    """
    Determine the density of air.

    Arguments:
        p : float -- pressure, [Pa]
        T : float -- temperature, [K]
        R : float -- specific gas constant, [J/(K*kg)]

    Returns:
        rho : float -- density, [kg/m^3]
    """

    # Determine the density
    rho = p / (R * T)

    return rho


def kinematic_viscosity(p, T, R, S, beta_s):
    """
    Determine the kinematic viscosity of air.

    Arguments:
        p : float -- pressure, [Pa]
        T : float -- temperature, [K]
        R : float -- specific gas constant, [J/(K*kg)]
        S : float -- sutherland's empirical coefficient, [K]
        beta_s : float -- sutherland's empirical coefficient, [kg/(m*s*K^0.5)]

    Returns:
        nu : float -- kinematic viscosity, [m^2/s]
    """

    # Determine the dynamic viscosity
    mu = beta_s * np.power(T, 3/2) / (T + S)

    # Determine the density
    rho = density(p, T, R)

    # Determine the kinematic viscosity
    nu = mu / rho

    return nu


def speed_of_sound(T, R, gamma):
    """
    Determine the speed of sound in air.

    Arguments:
        T : float -- temperature, [K]
        R : float -- specific gas constant, [J/(kg*K)]
        gamma : float -- specific heat ratio, [-]

    Returns:
        c : float -- speed of sound, [m/s]
    """

    # Determine the speed of sound
    c = np.sqrt(gamma * R * T)

    return c

if __name__ == "__main__":

    # Show the density of air
    print("Rho: " + str(density(101325.0, 288.15, 287.1)) + " [kg/m^3]")

    # Show the kinematic viscosity of air
    print("Nu: " + str(kinematic_viscosity(101325.0, 288.15, 287.1, 110.4, 1.458E-6)) + " [m^2/s]")

    # Show the speed of sound in air
    print("C: " + str(speed_of_sound(288.15, 287.1, 1.4)) + " [m/s]")
