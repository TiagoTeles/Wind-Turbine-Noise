"""
Author:   T. Moreira da Fonte Fonseca Teles
Email:    tmoreiradafont@tudelft.nl
Date:     2025-11-11
License:  GNU GPL 3.0

Store the salinity data.

Classes:
    Salinity 

Functions:
    None

Exceptions:
    None
"""

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import scipy as sp

from source.settings import LAT, LON, SALINITY_PATH


class Salinity:
    """
    A class to store the salinity data.

    Methods:
        __init__ -- initialise the Salinity class
        get -- get the salinity at a given latitude, longitude, and altitude
        show -- show the salinity profile at a given latitude and longitude

    Attributes:
        path : str -- path to the salinity file
        latitude : np.ndarray -- latitude, [rad]
        longitude : np.ndarray -- longitude, [rad]
        altitude : np.ndarray -- altitude, [m]
        salinity : np.ndarray -- salinity, [-]
        interpolator : RegularGridInterpolator -- 3D interpolator
    """

    def __init__(self, path):
        """
        Initialise the Salinity class.

        Parameters:
            path : str -- path to the salinity file

        Returns:
            None
        """

        self.path = path

        # Read the file
        data = pd.read_csv(self.path)

        # Determine the number latitudes and longitudes
        n_latitude = len(np.unique(data["latitude"]))
        n_longitude = len(np.unique(data["longitude"]))
        n_altitude = len(np.unique(data["altitude"]))

        # Reshape the data from a 1D array into a 2D array
        self.latitude = np.reshape(data["latitude"], (n_latitude, n_longitude, n_altitude))
        self.longitude = np.reshape(data["longitude"], (n_latitude, n_longitude, n_altitude))
        self.altitude = np.reshape(data["altitude"], (n_latitude, n_longitude, n_altitude))
        self.salinity = np.reshape(data["salinity"], (n_latitude, n_longitude, n_altitude))

        # Convert the coordinates from [deg] to [rad]
        self.latitude = np.radians(self.latitude)
        self.longitude = np.radians(self.longitude)

        # Create a 3D interpolator
        latitude = self.latitude[:, 0, 0]
        longitude = self.longitude[0, :, 0]
        altitude = self.altitude[0, 0, :]

        self.interpolator = sp.interpolate.RegularGridInterpolator((latitude, longitude, altitude), self.salinity)

    def get(self, latitude, longitude, altitude):
        """
        Get the salinity at a given latitude, longitude, and altitude.

        Parameters:
            latitude : np.ndarray -- latitude, [rad]
            longitude : np.ndarray -- longitude, [rad]
            altitude : np.ndarray -- altitude, [m]

        Returns:
            salinity : np.ndarray -- salinity, [-]
        """

        # Determine the salinity
        salinity = self.interpolator((latitude, longitude, altitude))

        return salinity

    def show(self, latitude, longitude, s_min, s_max):
        """
        Show the salinity profile at a given latitude and longitude.
        
        Parameters:
            latitude : float -- latitude, [rad]
            longitude : float -- longitude, [rad]
            s_min : float -- minimum salinity, [-]
            s_max : float -- maximum salinity, [-]
            
        Returns:
            None
        """

        # Determine the altitude array
        altitude = self.altitude[0, 0, :]

        # Determine the latitude and longitude arrays
        latitude = latitude * np.ones(altitude.shape)
        longitude = longitude * np.ones(altitude.shape)

        # Determine the salinity profile
        salinity = self.get(latitude, longitude, altitude)

        # Remove NaN values
        altitude = altitude[np.invert(np.isnan(salinity))]
        salinity = salinity[np.invert(np.isnan(salinity))]

        # Plot the salinity profile
        plt.plot(salinity, altitude)

        # Set the axis labels
        plt.xlabel("Salinity, [-]")
        plt.ylabel("Altitude, [m]")

        # Set the axis limits
        plt.xlim(s_min, s_max)
        plt.ylim(np.min(altitude), 0.0)

        # Use scientific notation on the x-axis
        plt.ticklabel_format(axis="x", style="scientific", scilimits=(-3, -3))

        # Show the plot
        plt.grid()
        plt.show()

if __name__ == "__main__":

    # Show the salinity profile
    salinity = Salinity(SALINITY_PATH)
    salinity.show(LAT, LON, 0.035, 0.036)
