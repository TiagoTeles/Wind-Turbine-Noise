"""
Author:   T. Moreira da Fonte Fonseca Teles
Email:    tmoreiradafont@tudelft.nl
Date:     2025-09-17
License:  GNU GPL 3.0

Store the seabed data.

Classes:
    Seabed 

Functions:
    None

Exceptions:
    None
"""

import matplotlib.patches as mpatches
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import scipy as sp

from source.ocean.bathymetry import Bathymetry
from source.settings import BATHYMETRY_PATH, SEABED_PATH


LABELS = {
    0: "Mud",
    1: "Sandy Mud",
    2: "Muddy Sand",
    3: "Sand",
    4: "Mixed Sediment",
    5: "Coarse Substrate",
    6: "Rock and Boulders",
}

class Seabed:
    """
    A class to store the seabed data.

    Methods:
        __init__ -- initialise the Seabed class
        get -- get the seabed type at a given latitude and longitude
        show -- show the seabed data

    Attributes:
        path : str -- path to the seabed file
        latitude : np.ndarray -- latitude, [rad]
        longitude : np.ndarray -- longitude, [rad]
        type : np.ndarray -- seabed type
        interpolator : RegularGridInterpolator -- 2D interpolator
    """

    def __init__(self, path):
        """
        Initialise the Seabed class.

        Parameters:
            path : str -- path to the seabed file

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
        self.type = np.reshape(data["type"], (n_latitude, n_longitude))

        # Convert the coordinates from [deg] to [rad]
        self.latitude = np.radians(self.latitude)
        self.longitude = np.radians(self.longitude)

        # Create a 2D interpolator
        latitude = self.latitude[:, 0]
        longitude = self.longitude[0, :]

        self.interpolator = sp.interpolate.RegularGridInterpolator((latitude, longitude), self.type, method="nearest")

    def get(self, latitude, longitude):
        """
        Get the seabed type at a given latitude and longitude.

        Parameters:
            latitude : np.array -- latitude, [rad]
            longitude : np.array -- longitude, [rad]

        Returns:
            type : np.array -- seabed type
        """

        # Determine the seabed type
        type = self.interpolator((latitude, longitude))

        return type

    def show(self, c_map, bathymetry):
        """
        Show the seabed data.

        Parameters:
            c_map : str -- colormap

        Returns:
            None
        """

        # Determine the image extent
        lat_min = np.degrees(np.min(self.latitude))
        lat_max = np.degrees(np.max(self.latitude))
        lon_min = np.degrees(np.min(self.longitude))
        lon_max = np.degrees(np.max(self.longitude))

        extent = [lon_min, lon_max, lat_min, lat_max]

        # Plot the seabed type
        v_min = 0
        v_max = plt.get_cmap(c_map).N - 1

        img = plt.imshow(self.type, cmap=c_map, vmin=v_min, vmax=v_max, origin="lower", extent=extent)

        # Set the legend
        patches = []

        for key, label in LABELS.items():
            color = img.cmap(img.norm(key))
            patches.append(mpatches.Patch(color=color, label=label))

        plt.legend(handles=patches)

        # Plot the coastline
        latitude = np.degrees(bathymetry.latitude)
        longitude = np.degrees(bathymetry.longitude)
        elevation = np.nan_to_num(bathymetry.elevation, nan=np.nextafter(0.0, 1.0))

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

    # Get the bathymetry
    bathymetry = Bathymetry(BATHYMETRY_PATH)

    # Show the seabed type
    seabed = Seabed(SEABED_PATH)
    seabed.show("tab10", bathymetry)
