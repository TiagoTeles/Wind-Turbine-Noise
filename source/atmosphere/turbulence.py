"""
Author:   T. Moreira da Fonte Fonseca Teles
Email:    tmoreiradafont@tudelft.nl
Date:     2025-11-12
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

from source.atmosphere.ISA import density, kinematic_viscosity
from source.constants import G, P, ALPHA_CH, ALPHA_M, KAPPA
from source.settings import P_0, T_0


def turbulence_intensity(z, z_0):
    """
    Determine the turbulence intensity.

    Parameters:
        z : np.ndarray -- height, [m]
        z_0 : np.ndarray -- surface roughness length, [m]

    Returns:
        I : np.ndarray -- turbulence intensity, [-]
    """

    # Determine the turbulence intensity
    gamma = 0.24 + 0.096 * np.log10(z_0) + 0.016 * np.square(np.log10(z_0))
    I = gamma * np.log10(30.0/z_0) / np.log10(z/z_0)

    return I


def turbulence_length_scale(z, z_0):
    """
    Determine the turbulence length scale.

    Parameters:
        z : np.ndarray -- height, [m]
        z_0 : np.ndarray -- surface roughness length, [m]

    Returns:
        L : np.ndarray -- turbulence length scale, [m]
    """

    # Determine the turbulence length scale
    L = 25.0 * np.pow(z, 0.35) * np.pow(z_0, -0.063)

    return L


def surface_roughness_length(u_n, z, nu, output_all=False):
    """
    Determine the surface roughness length.

    Parameters:
        u_n : np.ndarray -- neutral wind speed, [m/s]
        z : np.ndarray -- height, [m]
        nu : np.ndarray -- kinematic viscosity, [m^2/s]
        output_all : bool -- output all contributions?

    Returns:
        z_0 : np.ndarray -- surface roughness length, [m]
    """

    # Determine R and A
    R = z / (ALPHA_M * nu) * (KAPPA * u_n)
    A = ALPHA_CH / (G * z) * np.square(KAPPA * u_n)

    # Determine b_n
    b_n_nu = -1.47 + 0.93 * np.log(R)
    b_n_alpha = 2.65 - 1.44 * np.log(A) - 0.015 * np.square(np.log(A))
    b_n_fit = np.pow(np.pow(b_n_nu, P) + np.pow(b_n_alpha, P), 1.0 / P)

    # Determine the surface roughness length
    z_0_nu = z / (np.exp(b_n_nu) - 1.0)
    z_0_alpha = z / (np.exp(b_n_alpha) - 1.0)
    z_0 = z / (np.exp(b_n_fit) - 1.0)

    if output_all:
        return z_0_nu, z_0_alpha, z_0
    else:
        return z_0

if __name__ == "__main__":

    # Show the surface roughness length
    rho = density(P_0, T_0)
    nu = kinematic_viscosity(T_0, rho)
    z_0 = surface_roughness_length(11.0, 170.0, nu)

    print(f"Surface Roughness Length: {z_0} [m]")

    # Show the turbulence intensity
    z = np.linspace(1.0E-3, 200.0, 1000)
    I = turbulence_intensity(z, z_0)

    plt.plot(I, z)
    plt.xlabel("Turbulence Intensity, [-]")
    plt.ylabel("Height, [m]")
    plt.xlim(0.0, 0.2)
    plt.ylim(0.0, 200.0)
    plt.grid()
    plt.show()

    # Show the turbulence length scale
    z = np.linspace(0.0, 200.0, 1000)
    L = turbulence_length_scale(z, z_0)

    plt.plot(L, z)
    plt.xlabel("Turbulence Length Scale, [m]")
    plt.ylabel("Height, [m]")
    plt.xlim(0.0, 300.0)
    plt.ylim(0.0, 200.0)
    plt.grid()
    plt.show()

    # Show the surface roughness length
    u_n = np.linspace(1.0E-3, 12.0, 1000)
    z_0_nu, z_0_alpha, z_0 = surface_roughness_length(u_n, 170.0, nu, True)

    plt.semilogy(u_n, z_0_nu, label="Light Wind")
    plt.semilogy(u_n, z_0_alpha, label="Strong Wind")
    plt.semilogy(u_n, z_0, label="Empirical Fit")
    plt.xlabel("Neutral Wind Speed, [m/s]")
    plt.ylabel("Surface Roughness Length, [m]")
    plt.xlim(0.0, 12.0)
    plt.ylim(1.0E-5, 1.0E-3)
    plt.grid(which="both")
    plt.legend()
    plt.show()
