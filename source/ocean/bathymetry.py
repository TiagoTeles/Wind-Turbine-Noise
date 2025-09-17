"""
Author:   T. Moreira da Fonte Fonseca Teles
Email:    tmoreiradafont@tudelft.nl
Date:     2025-09-17
License:  GNU GPL 3.0

Store the bathymetry data.

Classes:
    Bathymetry 

Functions:
    None

Exceptions:
    None
"""

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import scipy as sp

from source.settings import BATHYMETRY_PATH


class Bathymetry:
    """
    A class to store the bathymetry data.

    Methods:
        __init__ -- initialise the Bathymetry class
        get -- get the elevation at a given latitude and longitude
        show -- show the bathymetry data

    Attributes:
        path : str -- path to the bathymetry file
        latitude : np.ndarray -- latitude, [rad]
        longitude : np.ndarray -- longitude, [rad]
        elevation : np.ndarray -- elevation, [m]
        interpolator : RegularGridInterpolator -- 2D interpolator
    """

    def __init__(self, path):
        """
        Initialise the Bathymetry class.

        Parameters:
            path : str -- path to the bathymetry file

        Returns:
            None
        """

        self.path = path

        # Read the file
        data = pd.read_csv(self.path)

        # Determine the number latitudes and longitudes
        n_latitude = len(np.unique(data["latitude"]))
        n_longitude = len(np.unique(data["longitude"]))

        # Reshape the data from a 1D array into a 2D array
        self.latitude = np.reshape(data["latitude"], (n_latitude, n_longitude))
        self.longitude = np.reshape(data["longitude"], (n_latitude, n_longitude))
        self.elevation = np.reshape(data["elevation"], (n_latitude, n_longitude))

        # Convert the coordinates from [deg] to [rad]
        self.latitude = np.radians(self.latitude)
        self.longitude = np.radians(self.longitude)

        # Create a 2D interpolator
        latitude = self.latitude[:, 0]
        longitude = self.longitude[0, :]

        self.interpolator = sp.interpolate.RegularGridInterpolator((latitude, longitude), self.elevation)

    def get(self, latitude, longitude):
        """
        Get the elevation at a given latitude and longitude.

        Parameters:
            latitude : np.array -- latitude, [rad]
            longitude : np.array -- longitude, [rad]

        Returns:
            elevation : np.array -- elevation, [m]
        """

        # Determine the elevation
        elevation = self.interpolator((latitude, longitude))

        return elevation

    def show(self, c_map, z_min, z_max):
        """
        Show the bathymetry data.

        Parameters:
            c_map : str -- colormap
            z_min : float -- minimum elevation
            z_max : float -- maximum elevation

        Returns:
            None
        """

        # Determine the image extent
        lat_min = np.degrees(np.min(self.latitude))
        lat_max = np.degrees(np.max(self.latitude))
        lon_min = np.degrees(np.min(self.longitude))
        lon_max = np.degrees(np.max(self.longitude))

        extent = [lon_min, lon_max, lat_min, lat_max]

        # Plot the elevation
        plt.imshow(self.elevation, cmap=c_map, vmin=z_min, vmax=z_max, origin="lower", extent=extent)
        plt.colorbar(label="Elevation, [m]", extend="min")

        # Plot the coastline
        latitude = np.degrees(self.latitude)
        longitude = np.degrees(self.longitude)
        elevation = np.nan_to_num(self.elevation, nan=np.nextafter(0.0, 1.0))

        plt.contour(longitude, latitude, elevation, levels=[0.0], colors="black")

        # Set the axis labels
        plt.xlabel(r"Longitude, [$^\circ$]")
        plt.ylabel(r"Latitude, [$^\circ$]")

        # Set the axis limits
        plt.xlim(lon_min, lon_max)
        plt.ylim(lat_min, lat_max)

        # Show the plot
        plt.show()

if __name__ == "__main__":

    # Show the bathymetry
    bathymetry = Bathymetry(BATHYMETRY_PATH)
    bathymetry.show("inferno", -200.0, 0.0)
