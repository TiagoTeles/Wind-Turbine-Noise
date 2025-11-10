"""
Author:   T. Moreira da Fonte Fonseca Teles
Email:    tmoreiradafont@tudelft.nl
Date:     2025-11-10
License:  GNU GPL 3.0

Determine the sound speed profile.

Classes:
    None 

Functions:
    pressure
    sound_speed_profile

Exceptions:
    None
"""

import matplotlib.pyplot as plt
import numpy as np

from source.constants import G
from source.ocean.salinity import Salinity
from source.ocean.temperature import Temperature
from source.settings import P_0, LAT, LON, SALINITY_PATH, TEMPERATURE_PATH


def pressure(z, p_0, rho):
    """
    Determine the ocean pressure.

    Parameters:
        z : np.ndarray -- altitude, [m]
        p_0 : float -- pressure, [Pa]
        rho : np.ndarray -- density, [kg/m^3]

    Returns:
        p : np.ndarray -- pressure, [Pa]
    """

    # Check for positive altitudes
    if np.any(z > 0.0):
        print("Positive altitude detected! z > 0.0 [m].")

    # Determine the pressure
    p = p_0 - rho * G * z

    return p


def sound_speed_profile(T, S, P):
    """
    Determine the sound speed profile.

    Parameters:
        T : np.ndarray -- temperature, [K]
        S : np.ndarray -- salinity, [-]
        P : np.ndarray -- pressure, [Pa]

    Returns:
        c : np.ndarray -- sound speed, [m/s]
    """

    # Check for low temperatures
    if np.any(T < 273.15):
        print("Low ocean temperature detected! T < 273.15 [K].")

    # Check for high temperatures
    if np.any(T > 313.15):
        print("High ocean temperature detected! T > 313.15 [K].")

    # Check for high salinities
    if np.any(S > 0.04):
        print("High ocean salinity detected! S > 0.04 [-].")

    # Check for high pressures
    if np.any(P > 1.0E8):
        print("High ocean pressure detected! P > 1.00E8 [Pa].")

    # Convert the temperature from [K] to [°C]
    T = T - 273.15

    # Convert the salinity from [-] to [ppt]
    S = S * 1E3

    # Convert the pressure from [Pa] to [bar]
    P = P / 1E5

    # Define the UNESCO model coefficients
    c_coeffs = np.array([[  1402.388,    5.03830, -5.81090E-2,   3.3432E-4, -1.47797E-6, 3.1419E-9],
                         [  0.153563,  6.8999E-4,  -8.1829E-6,   1.3632E-7, -6.1260E-10,       0.0],
                         [ 3.1260E-5, -1.7111E-6,   2.5986E-8, -2.5353E-10,  1.0415E-12,       0.0],
                         [-9.7729E-9, 3.8513E-10, -2.3654E-12,         0.0,         0.0,       0.0]])

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
    c = C_w + A * S + B * np.pow(S, 3/2) + D * np.pow(S, 2)

    return c

if __name__ == "__main__":

    # Read the study domain data
    salinity = Salinity(SALINITY_PATH)
    temperature = Temperature(TEMPERATURE_PATH)

    # Determine the altitude, latitude, longitude 
    altitude = salinity.altitude[0, 0, :]
    latitude = LAT * np.ones(altitude.shape)
    longitude = LON * np.ones(altitude.shape)

    # Determine the salinity and temperature
    S = salinity.get(latitude, longitude, altitude)
    T = temperature.get(latitude, longitude, altitude)

    # Determine the pressure
    P = pressure(altitude, P_0, 1025.0)

    # Determine the sound speed profile
    c = sound_speed_profile(T, S, P)

    # Remove NaN values
    altitude = altitude[np.invert(np.isnan(c))]
    c = c[np.invert(np.isnan(c))]

    # Show the sound speed profile
    plt.plot(c, altitude)
    plt.xlabel("Speed of Sound, [m/s]")
    plt.ylabel("Altitude, [m]")
    plt.xlim(1502.0, 1514.0)
    plt.ylim(np.min(altitude), 0.0)
    plt.grid()
    plt.show()
