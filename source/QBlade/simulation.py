"""
Author:   T. Moreira da Fonte Fonseca Teles
Email:    tmoreiradafont@tudelft.nl
Date:     2025-09-22
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

import pandas as pd

from source.QBlade.turbine import Turbine


class Simulation:
    """
    A class to store the simulation data.

    Methods:
        __init__ -- initialise the Simulation class
        read_results -- read the simulation results

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

        # Add the Turbine object
        turbine_path = os.path.normpath(os.path.join(os.path.dirname(self.path), lines[16].split()[0]))
        self.turbine = Turbine(turbine_path)

        # Close the file
        f.close()

        # Initialise the results
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

        # Read the Simulation results
        self.time = self.results["Time~[s]"].to_numpy()
        self.inflow_velocity = self.results["Abs~Inflow~Vel.~at~Hub~[m/s]"].to_numpy()

        # Read the Turbine results
        self.turbine.read_results(self.results)
