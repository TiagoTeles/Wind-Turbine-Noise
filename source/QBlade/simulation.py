"""
Author:   T. Moreira da Fonte Fonseca Teles
Email:    tmoreiradafont@tudelft.nl
Date:     2025-11-11
License:  GNU GPL 3.0

Store the simulation data.

Classes:
    Simulation

Functions:
    None

Exceptions:
    None
"""

import os
import sys

import numpy as np
import pandas as pd

from source.QBlade.turbine import Turbine


class Simulation:
    """
    A class to store the simulation data.

    Methods:
        __init__ -- initialise the Simulation class
        read_results -- read the simulation results
        get_results -- get the simulation results

    Attributes:
        path : str -- path to the .sim file
        turbine: Turbine -- turbine object
        results : pd.DataFrame -- simulation results
        time : np.ndarray -- time, [s]
        inflow_velocity : np.ndarray -- inflow velocity, [m/s]
    """

    def __init__(self, path):
        """
        Initialise the Simulation class.

        Parameters:
            path : str -- path to the .sim file

        Returns:
            None
        """

        self.path = path

        # Read the file
        f = open(self.path, "r", encoding="utf-8")
        lines = f.readlines()

        # Determine the Turbine path
        turbine_directory = os.path.dirname(self.path)
        turbine_name = lines[16].split()[0]
        turbine_path = os.path.normpath(os.path.join(turbine_directory, turbine_name))

        # Add the Turbine object
        self.turbine = Turbine(turbine_path)

        # Close the file
        f.close()

        # Initialise the results attributes
        self.results = None
        self.time = None
        self.inflow_velocity = None

    def read_results(self, path):
        """
        Read the simulation results.

        Parameters:
            path : str -- path to the results file

        Returns:
            None
        """

        # Read the results
        self.results = pd.read_csv(path, delimiter=r"\s+", skiprows=2)

        # Determine the start and end indices of the last revolution
        index_0 = np.where(np.diff(self.results["Azimuthal~Angle~BLD_1~[deg]"]) < 0.0)[0][-2] + 1
        index_1 = np.where(np.diff(self.results["Azimuthal~Angle~BLD_1~[deg]"]) < 0.0)[0][-1] + 1

        # Keep the last full revolution
        self.results = self.results.iloc[index_0:index_1]

        # Read the Simulation results
        self.time = self.results["Time~[s]"].to_numpy()
        self.inflow_velocity = self.results["Abs~Inflow~Vel.~at~Hub~[m/s]"].to_numpy()

        # Read the Turbine results
        self.turbine.read_results(self.results)

    def get_results(self, key, azimuth):
        """
        Get the simulation results.

        Parameters:
            key : str -- member name
            azimuth : float -- azimuth angle, [rad]
        """

        # Check if the key is valid
        if key in ["time", "inflow_velocity"]:

            # Interpolate the results
            value = np.interp(azimuth, self.turbine.azimuth, getattr(self, key), period=2*np.pi)

        else:
            print("Invalid key!")
            sys.exit(1)

        return value
