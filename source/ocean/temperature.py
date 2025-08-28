"""
Author:   T. Moreira da Fonte Fonseca Teles
Email:    tmoreiradafont@tudelft.nl
Date:     2025-08-25
License:  GNU GPL 3.0

Store the temperature data.

Classes:
    Temperature 

Functions:
    None

Exceptions:
    None
"""

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import scipy as sp


class Temperature:
    """
    A class to store the temperature data.

    Methods:
        __init__ -- initialise the Temperature class
        get -- get the temperature at a given latitude, longitude, and altitude
        show -- show the temperature profile at a given latitude and longitude

    Attributes:
        altitude : np.ndarray -- altitude, [m]
        interpolator : RegularGridInterpolator -- 3D interpolator
        latitude : np.ndarray -- latitude, [rad]
        longitude : np.ndarray -- longitude, [rad]
        path : str -- path to the temperature file
    """

    def __init__(self, path):
        """
        Initialise the Temperature class.

        Arguments:
            path : str -- path to the temperature file

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
        self.temperature = np.reshape(data["temperature"], (n_latitude, n_longitude, n_altitude))

        # Convert the coordinates from [deg] to [rad]
        self.latitude = np.radians(self.latitude)
        self.longitude = np.radians(self.longitude)

        # Create a 3D interpolator
        latitude = self.latitude[:, 0, 0]
        longitude = self.longitude[0, :, 0]
        altitude = self.altitude[0, 0, :]

        self.interpolator = sp.interpolate.RegularGridInterpolator((latitude, longitude, altitude), self.temperature)

    def get(self, latitude, longitude, altitude):
        """
        Get the temperature at a given latitude, longitude, and altitude.

        Arguments:
            latitude : np.array -- latitude, [rad]
            longitude : np.array -- longitude, [rad]
            altitude : np.array -- altitude, [m]

        Returns:
            temperature : np.array -- temperature, [-]
        """

        coordinate = np.stack((latitude, longitude, altitude), axis=1)

        temperature = self.interpolator(coordinate)

        return temperature

    def show(self, latitude, longitude, t_min, t_max):
        """
        Show the temperature profile at a given latitude and longitude.
        
        Arguments:
            latitude : float -- latitude, [rad]
            longitude : float -- longitude, [rad]
            t_min : float -- minimum temperature, [°C]
            t_max : float -- maximum temperature, [°C]

        Returns:
            None
        """

        # Determine the altitude array
        altitude = self.altitude[0, 0, :]

        # Determine the latitude and longitude arrays
        latitude = latitude * np.ones(len(altitude))
        longitude = longitude * np.ones(len(altitude))

        # Determine the temperature profile
        temperature = self.get(latitude, longitude, altitude)

        # Remove NaN values
        altitude = altitude[np.invert(np.isnan(temperature))]
        temperature = temperature[np.invert(np.isnan(temperature))]

        # Plot the temperature profile
        plt.plot(temperature - 273.15, altitude)

        # Set the axis labels
        plt.xlabel(r"Temperature, [$^\circ C$]")
        plt.ylabel("Altitude, [m]")

        # Set the axis limits
        plt.xlim(t_min, t_max)
        plt.ylim(np.min(altitude), 0)

        # Show the plot
        plt.grid()
        plt.show()

if __name__ == "__main__":

    # Show the temperature profile
    temperature = Temperature("data\\environments\\windfloat_atlantic\\temperature.csv")
    temperature.show(np.radians(41.6865), np.radians(-9.0574), 13, 17)
