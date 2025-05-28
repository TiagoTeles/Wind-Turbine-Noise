"""
Author:   T. Moreira da Fonte Fonseca Teles
Email:    tmoreiradafont@tudelft.nl
Date:     2025-05-05
License:  GNU GPL 3.0

Miscellaneous functions.

Classes:
    None

Functions:
    octave
    turbulence_intensity
    turbulence_length_scale
    surface_roughness_length
    gauge_pressure
    sound_speed_profile

Exceptions:
    None
"""

import numpy as np


def octave(f_min, f_max, f_ref, base_10):
    """
    Determine the center, lower, and upper frequencies of the 1/3 octave frequency bands.

    Arguments:
        f_min : float -- Minimum frequency, [Hz]
        f_max : float -- Maximum frequency, [Hz]
        f_ref : float -- Reference frequency, [Hz]
        base_10 : bool -- Use base 10?

    Returns:
        f_center : np.array -- Center frequencies, [Hz]
        f_lower : np.array -- Lower frequencies, [Hz]
        f_upper : np.array -- Upper frequencies, [Hz]
    """

    if base_10:
        # Determine the smallest and largest index
        min_index = np.floor(10 * np.log10(f_min / f_ref))
        max_index = np.ceil(10 * np.log10(f_max / f_ref))

        # Determine the center, lower and upper frequencies
        f_center = f_ref * np.pow(10, np.arange(min_index, max_index+1) / 10)
        f_lower = f_center / np.pow(10, 1/20)
        f_upper = f_center * np.pow(10, 1/20)

    else:
        # Determine the smallest and largest index
        min_index = np.floor(3 * np.log2(f_min / f_ref))
        max_index = np.ceil(3 * np.log2(f_max / f_ref))

        # Determine the center, lower and upper frequencies
        f_center = f_ref * np.pow(2, np.arange(min_index, max_index+1) / 3)
        f_lower = f_center / np.pow(2, 1/6)
        f_upper = f_center * np.pow(2, 1/6)

    return f_center, f_lower, f_upper


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

    # LOW SURFACE ROUGHNESS DETECTED!

    return L


def surface_roughness_length(z_ref, U_ref, g, kappa, nu):
    """
    Determine the surface roughness length.

    Parameters:
        z_ref : np.array -- reference height, [m]
        U_ref : np.array -- reference velocity, [m/s]
        g : np.array -- gravitational acceleration, [m/s^2]
        kappa : np.array -- Von Karman constant, [-]
        nu : np.array -- kinematic viscosity, [m^2/s]

    Returns:
        z_0 : np.array -- surface roughness length, [m]
    """

    # Determine R and A
    R = z_ref / (0.11 * nu) * (kappa * U_ref)
    A = 0.018 / (g * z_ref) * np.square(kappa * U_ref)

    # Determine b_n
    b_n_nu = -1.47 + 0.93 * np.log(R)
    b_n_alpha = 2.65 - 1.44 * np.log(A) - 0.015 * np.square(np.log(A))
    b_n = np.pow(np.pow(b_n_nu, -12) + np.pow(b_n_alpha, -12), -1/12)

    # Determine the roughness length
    z_0 = z_ref / (np.exp(b_n) - 1)

    return z_0


def gauge_pressure(Z, phi):
    """
    Determine the ocean pressure.

    Parameters:
        Z : np.array -- depth, [m]
        phi : np.array -- latitude, [rad]

    Returns:
        P : np.array -- gauge pressure, [Pa]
    """

    # Check the high latitudes
    if np.degrees(phi) > 60:
        print("High latitude detected! phi > 60 [°].")

    # Check for low latitudes
    if np.degrees(phi) < -40:
        print("Low latitude detected! phi < -40 [°].")

    # Determine the pressure for a standard ocean
    h_45 = 1.00818E-2 * Z + 2.465E-8 * np.square(Z) \
         - 1.25E-13 * np.power(Z, 3) + 2.8E-19 * np.power(Z, 4)

    g = 9.7803 * (1 + 5.3E-3 * np.square(np.sin(phi)))

    k = (g - 2E-5 * Z) / (9.80612 - 2E-5 * Z)

    h = h_45 * k

    # Determine the pressure correction
    h_0_Z = 1E-2 * Z / (Z + 100) + 6.2E-6 * Z

    # Determine the pressure
    P = h - h_0_Z

    # Convert the pressure from [MPa] to [Pa]
    P *= 1E6

    return P


def sound_speed_profile(T, S, P):
    """
    Determine the sound speed profile.

    Parameters:
        T : np.array -- ocean temperature, [K]
        S : np.array -- ocean salinity, [-]
        P : np.array -- ocean pressure, [Pa]

    Returns:
        c : np.array -- sound speed, [m/s]

    """

    # Check for low ocean temperatures
    if T < 273.15:
        print("Low ocean temperature detected! T < 273.15 [K].")

    # Check for high ocean temperatures
    if T > 313.15:
        print("High ocean temperature detected! T > 313.15 [K].")

    # Check for high ocean salinities
    if S > 0.04:
        print("High ocean salinity detected! S > 0.04 [-].")

    # Check for high ocean pressures
    if P > 1E8:
        print("High ocean pressure detected! P > 1.00E8 [Pa].")

    # Convert the temperature from [K] to [°C]
    T -= 273.15

    # Convert the salinity from [-] to [ppt]
    S *= 1E3

    # Convert the pressure from [Pa] to [bar]
    P /= 1E5

    # Define the UNESCO model coefficients
    c_coeffs = np.array([[  1402.388,    5.03830, -5.81090E-2,   3.3432E-4, -1.47797E-6, 3.1419E-9],
                         [  0.153563,  6.8999E-4,  -8.1829E-6,   1.3632E-7, -6.1260E-10,       0.0],
                         [ 3.1260E-5, -1.7111E-6,   2.5986E-8, -2.5353E-10,  1.0415E-12,       0.0],
                         [-9.7729E-9,  3.8513E-1, -2.3654E-12,         0.0,         0.0,       0.0]])

    a_coeffs = np.array([[     1.389,  -1.262E-2,    7.166E-5,  2.008E-6,    -3.21E-8],
                         [ 9.4742E-5, -1.2583E-5,  -6.4928E-8, 1.0515E-8, -2.0142E-10],
                         [-3.9064E-7,  9.1061E-9, -1.6009E-10, 7.994E-12,         0.0],
                         [ 1.100E-10,  6.651E-12,  -3.391E-13,       0.0,         0.0]])

    b_coeffs = np.array([[-1.922E-2,  -4.42E-5],
                         [7.3637E-5, 1.7950E-7]])

    d_coeffs = np.array([[  1.727E-3],
                         [-7.9836E-6]])

    # Determine C_w, A, B, and D
    C_w = np.polynomial.polynomial.polyval2d(P, T, c_coeffs)
    A = np.polynomial.polynomial.polyval2d(P, T, a_coeffs)
    B = np.polynomial.polynomial.polyval2d(P, T, b_coeffs)
    D = np.polynomial.polynomial.polyval2d(P, T, d_coeffs)

    # Determine the sound speed profile
    c = C_w + A * S + B * np.power(S, 3/2) + D * np.power(S, 2)

    return c
