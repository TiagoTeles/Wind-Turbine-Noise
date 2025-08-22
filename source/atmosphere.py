"""
Author:   T. Moreira da Fonte Fonseca Teles
Email:    tmoreiradafont@tudelft.nl
Date:     2025-08-22
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
        p : float -- Pressure, [Pa]
        T : float -- Temperature, [K]
        R : float -- Specific gas constant, [J/(K*kg)]

    Returns:
        rho : float -- Density, [kg/m^3]
    """

    # Determine the density
    rho = p / (R * T)

    return rho


def kinematic_viscosity(p, T, R, S, beta_s):
    """
    Determine the kinematic viscosity of air.

    Arguments:
        p : float -- Pressure, [Pa]
        T : float -- Temperature, [K]
        R : float -- Specific gas constant, [J/(K*kg)]
        S : float -- Sutherland's empirical coefficient, [K]
        beta_s : float -- Sutherland's empirical coefficient, [kg/(m*s*K^0.5)]

    Returns:
        nu : float -- Kinematic viscosity, [m^2/s]
    """

    # Determine the dynamic viscosity
    mu = beta_s * np.power(T, 1.5) / (T + S)

    # Determine the density
    rho = density(p, T, R)

    # Determine the kinematic viscosity
    nu = mu / rho

    return nu


def speed_of_sound(T, R, gamma):
    """
    Determine the speed of sound in air.

    Arguments:
        T : float -- Temperature, [K]
        R : float -- Specific gas constant, [J/(kg*K)]
        gamma : float -- Specific heat ratio, [-]

    Returns:
        c : float -- Speed of sound, [m/s]
    """

    # Determine the speed of sound
    c = np.sqrt(gamma * R * T)

    return c
