"""
Author:   T. Moreira da Fonte Fonseca Teles
Email:    tmoreiradafont@tudelft.nl
Date:     2025-09-22
License:  GNU GPL 3.0

Store the simulation results.

Classes:
    Results

Functions:
    None

Exceptions:
    None
"""

import sys

import numpy as np
import pandas as pd


class Results:
    """
    A class to store the simulation results.

    Methods:
        __init__ -- initialise the Results class
        get -- get the results at a given azimuth angle

    Attributes:
        path : str -- path to the .txt file
        time : np.ndarray -- time, [s]
        velocity_hub : np.ndarray -- wind speed at the hub, [m/s]
        yaw : np.ndarray -- yaw angle, [rad]
        pitch : np.ndarray -- pitch angle, [rad]
        azimuth : np.ndarray -- azimuth angle, [rad]
        radius : np.ndarray -- spanwise position, [m]
        velocity : np.ndarray -- velocity, [m/s]
        aoa : np.ndarray -- angle of attack, [rad]
    """

    def __init__(self, path):
        """
        Initialise the Results class.

        Parameters:
            path : str -- path to the .txt file

        Returns:
            None
        """

        self.path = path

        # Read the simulation results
        results = pd.read_csv(self.path, delimiter=r"\s+", skiprows=2)

        # Determine the number of blades
        n_blades = 0

        for column in results.columns.values:
            if "Azimuthal~Angle~BLD_" in column:
                n_blades += 1

        # Determine the number of panels
        n_panels = 0

        for column in results.columns.values:
            if "Radius~BLD_1~PAN_" in column:
                n_panels += 1

        # Determine the turbine properties
        self.time = results["Time~[s]"].to_numpy()
        self.velocity_hub = results["Abs~Inflow~Vel.~at~Hub~[m/s]"].to_numpy()
        self.yaw = results["Yaw~Angle~[deg]"].to_numpy()

        # Determine the blade properties
        self.pitch = np.empty((len(results), n_blades))
        self.azimuth = np.empty((len(results), n_blades))

        for i in range(n_blades):
            self.pitch[:, i] = results[f"Pitch~Angle~BLD_{i+1}~[deg]"]
            self.azimuth[:, i] = results[f"Azimuthal~Angle~BLD_{i+1}~[deg]"]

        # Determine the panel properties
        self.radius = np.empty((len(results), n_blades, n_panels))
        self.velocity = np.empty((len(results), n_blades, n_panels))
        self.aoa = np.empty((len(results), n_blades, n_panels))

        for i in range(n_blades):
            for j in range(n_panels):
                self.radius[:, i, j] = results[f"Radius~BLD_{i+1}~PAN_{j}~[m]"]
                self.velocity[:, i, j] = results[f"Total~Velocity~BLD_{i+1}~PAN_{j}~[m/s]"]
                self.aoa[:, i, j] = results[f"Angle~of~Attack~at~0.25c~BLD_{i+1}~PAN_{j}~[deg]"]

        # Convert the angles from [deg] to [rad]
        self.yaw = np.radians(self.yaw)
        self.pitch = np.radians(self.pitch)
        self.azimuth = np.radians(self.azimuth)
        self.aoa = np.radians(self.aoa)

    def get(self, key, azimuth, blade, radius):
        """
        Get the results at a given time.

        Parameters:
            key : str -- property key
            azimuth : int -- azimuth angle, [rad]
            blade : int -- blade index
            radius : np.ndarray -- spanwise position, [m]
        
        Returns:
            value : np.ndarray -- property value
        """

        # Determine the time index
        offset = (self.azimuth[:, blade] - azimuth + np.pi) % (2 * np.pi) - np.pi
        index = np.where(np.diff(np.sign(offset)) > 0.0)[0][-1]

        # Determine the value
        if key == "time":
            value = self.time[index]

        elif key == "velocity_hub":
            value = self.velocity_hub[index]

        elif key == "yaw":
            value = self.yaw[index]

        elif key == "pitch":
            value = self.pitch[index, blade]

        elif key == "azimuth":
            value = self.azimuth[index, blade]

        elif key == "velocity":
            value = np.interp(radius, self.radius[index, blade, :], self.velocity[index, blade, :])

        elif key == "aoa":
            value = np.interp(radius, self.radius[index, blade, :], self.aoa[index, blade, :])

        else:
            print("Key not recognised!")
            sys.exit(1)

        return value
