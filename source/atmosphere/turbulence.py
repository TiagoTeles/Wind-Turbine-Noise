"""
Author:   T. Moreira da Fonte Fonseca Teles
Email:    tmoreiradafont@tudelft.nl
Date:     2025-08-28
License:  GNU GPL 3.0

Determine the turbulence characteristics.

Classes:
    None

Functions:
    turbulence_intensity
    turbulence_length_scale
    surface_roughness_length

Exceptions:
    None
"""

import matplotlib.pyplot as plt
import numpy as np

from source.atmosphere.international_standard_atmosphere import kinematic_viscosity
from source.constants import G, P, ALPHA_CH, ALPHA_M, KAPPA
from source.settings import P_0, T_0


def turbulence_intensity(z, z_0):
    """
    Determine the turbulence intensity.

    Parameters:
        z : np.array -- height above the terrain, [m]
        z_0 :  np.array -- surface roughness length, [m]

    Returns:
        I :  np.array -- turbulence intensity, [-]
    """

    gamma = 0.24 + 0.096 * np.log10(z_0) + 0.016 * np.square(np.log10(z_0))
    I = gamma * np.log10(30.0/z_0) / np.log10(z/z_0)

    return I


def turbulence_length_scale(z, z_0):
    """
    Determine the turbulence length scale.

    Parameters:
        z : np.array -- height above the terrain, [m]
        z_0 :  np.array -- surface roughness length, [m]

    Returns:
        L :  np.array -- turbulence length scale, [m]
    """

    L = 25.0 * np.pow(z, 0.35) * np.pow(z_0, -0.063)

    return L


def surface_roughness_length(z, u_n, nu, output_individual=False):
    """
    Determine the surface roughness length.

    Parameters:
        z : np.array -- height, [m]
        u_n : np.array -- neutral velocity, [m/s]
        nu : np.array -- kinematic viscosity, [m^2/s]
        output_individual : bool -- output individual contributions?

    Returns:
        z_0 : np.array -- surface roughness length, [m]
    """

    # Determine R and A
    R = z / (ALPHA_M * nu) * (KAPPA * u_n)
    A = ALPHA_CH / (G * z) * np.square(KAPPA * u_n)

    # Determine b_n
    b_n_nu = -1.47 + 0.93 * np.log(R)
    b_n_alpha = 2.65 - 1.44 * np.log(A) - 0.015 * np.square(np.log(A))
    b_n_fit = np.pow(np.pow(b_n_nu, P) + np.pow(b_n_alpha, P), 1.0/P)

    # Determine the surface roughness length
    z_0 = z / (np.exp(b_n_fit) - 1.0)

    # Determine the individual contributions
    z_0_nu = z / (np.exp(b_n_nu) - 1.0)
    z_0_alpha = z / (np.exp(b_n_alpha) - 1.0)

    # Return the surface roughness length
    if output_individual:
        return z_0, z_0_nu, z_0_alpha
    else:
        return z_0


if __name__ == "__main__":

    # Determine the kinematic viscosity of air
    nu = kinematic_viscosity(P_0, T_0)

    # Show the turbulence intensity
    z = np.linspace(1.0E-3, 200.0, 1000)
    z_0 = surface_roughness_length(170.0, 11.0, nu)
    I = turbulence_intensity(z, z_0)

    plt.plot(100.0 * I, z)
    plt.xlabel("Turbulence Intensity, [%]")
    plt.ylabel("Height, [m]")
    plt.xlim(0.0, 20.0)
    plt.ylim(0.0, 200.0)
    plt.grid()
    plt.show()

    # Show the turbulence length scale
    z = np.linspace(1.0E-3, 200.0, 1000)
    z_0 = surface_roughness_length(170.0, 11.0, nu)
    L = turbulence_length_scale(z, z_0)

    plt.plot(L, z)
    plt.xlabel("Turbulence Length Scale, [m]")
    plt.ylabel("Height, [m]")
    plt.xlim(0.0, 300.0)
    plt.ylim(0.0, 200.0)
    plt.grid()
    plt.show()

    # Show the surface roughness length
    U = np.linspace(1.0E-3, 12.0, 1000)
    z_0, z_0_nu, z_0_alpha = surface_roughness_length(170.0, U, nu, True)

    plt.semilogy(U, z_0_nu, label="Kinematic Viscosity")
    plt.semilogy(U, z_0_alpha, label="Charnock's Relation")
    plt.semilogy(U, z_0, label="Empirical Fit")
    plt.xlabel("Neutral Wind Speed, [m/s]")
    plt.ylabel("Surface Roughness Length, [m]")
    plt.xlim(0.0, 12.0)
    plt.ylim(1E-5, 1E-3)
    plt.grid(which="both")
    plt.legend()
    plt.show()
