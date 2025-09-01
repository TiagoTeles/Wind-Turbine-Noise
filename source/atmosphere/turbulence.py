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

from international_standard_atmosphere import kinematic_viscosity
from source.

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
    I = gamma * np.log10(30/z_0) / np.log10(z/z_0)

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

    L = 25 * np.pow(z, 0.35) * np.pow(z_0, -0.063)

    return L


def surface_roughness_length(z_ref, U_ref, nu, output_individual=False):
    """
    Determine the surface roughness length.

    Parameters:
        z_ref : np.array -- reference height, [m]
        U_ref : np.array -- reference velocity, [m/s]
        nu : np.array -- kinematic viscosity, [m^2/s]
        output_individual : bool -- output individual contributions?

    Returns:
        z_0 : np.array -- surface roughness length, [m]
    """

    # Determine R and A
    R = z_ref / (ALPHA_M * nu) * (KAPPA * U_ref)
    A = ALPHA_CH / (g * z_ref) * np.square(KAPPA * U_ref)

    # Determine b_n
    b_n_nu = -1.47 + 0.93 * np.log(R)
    b_n_alpha = 2.65 - 1.44 * np.log(A) - 0.015 * np.square(np.log(A))
    b_n = np.pow(np.pow(b_n_nu, p) + np.pow(b_n_alpha, p), 1.0/p)

    # Determine the surface roughness length
    z_0 = z_ref / (np.exp(b_n) - 1)

    # Determine the individual contributions
    z_0_nu = z_ref / (np.exp(b_n_nu) - 1)
    z_0_alpha = z_ref / (np.exp(b_n_alpha) - 1)

    # Return the surface roughness length
    if output_individual: 
        return z_0, z_0_nu, z_0_alpha
    else:
        return z_0


if __name__ == "__main__":

    # Determine the kinematic viscosity of air
    rho = density(101325.0, 288.15)
    nu = kinematic_viscosity(288.15, rho)

    # Show the turbulence intensity
    z = np.linspace(1.0E-3, 200.0, 1000)
    z_0 = surface_roughness_length(170.0, 11.0, nu)
    I = turbulence_intensity(z, z_0)

    plt.plot(100*I, z)
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
