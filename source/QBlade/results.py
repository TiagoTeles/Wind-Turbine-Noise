"""
Author:   T. Moreira da Fonte Fonseca Teles
Email:    tmoreiradafont@tudelft.nl
Date:     2025-09-17
License:  GNU GPL 3.0

Store the simulation results.

Classes:
    Results

Functions:
    None

Exceptions:
    None
"""

import os
import sys

import numpy as np
import pandas as pd


class Results:
    """
    A class to store the simulation results.

    Methods:
        __init__ -- initialise the Results class
        read -- read the results file

    Attributes:
        path : str -- path to the results file
        U_hub : ndarray -- wind speed at the hub, [m/s]
        phi : ndarray -- yaw angle, [rad]
        theta : ndarray -- pitch angle, [rad]
        psi : ndarray -- azimuth angle, [rad]
        r : ndarray -- spanwise position, [m]
        U : ndarray -- velocity, [m/s]
        alpha : ndarray -- angle of attack, [rad]
    """

    def __init__(self, path):
        """
        Initialise the Results class.

        Parameters:
            path : str -- path to the results file

        Returns:
            None
        """

        self.path = path

        # Check if the file exists
        if not os.path.isfile(path):
            print(f"No file found at {path}!")
            sys.exit(1)

        # Read the file
        self.read()

    def read(self):
        """
        Read the results file.

        Parameters:
            None

        Returns:
            None
        """

        # Read the simulation results
        data = pd.read_csv(self.path, delimiter=r"\s+", skiprows=2)

        # Determine the number of time steps
        n_steps = data.shape[0]

        # Determine the number of blades
        n_blades = 0

        for column in data.columns.values:
            if "Azimuthal~Angle~BLD_" in column:
                n_blades = n_blades + 1

        # Determine the number of panels
        n_panels = 0

        for column in data.columns.values:
            if "Radius~BLD_1~PAN_" in column:
                n_panels = n_panels + 1

        # Determine the turbine properties
        self.U_hub =  data["Abs~Inflow~Vel.~at~Hub~[m/s]"].to_numpy()
        self.phi = data["Yaw~Angle~[deg]"].to_numpy()

        # Determine the blade properties
        self.theta = np.empty((n_steps, n_blades))
        self.psi = np.empty((n_steps, n_blades))

        for i in range(n_blades):
            self.theta[:, i] = data[f"Pitch~Angle~BLD_{i+1}~[deg]"].to_numpy()
            self.psi[:, i] = data[f"Azimuthal~Angle~BLD_{i+1}~[deg]"].to_numpy()

        # Determine the panel properties
        self.r = np.empty((n_steps, n_blades, n_panels))
        self.U = np.empty((n_steps, n_blades, n_panels))
        self.alpha = np.empty((n_steps, n_blades, n_panels))

        for i in range(n_blades):
            for j in range(n_panels):
                self.r[:, i, j] = data[f"Radius~BLD_{i+1}~PAN_{j}~[m]"].to_numpy()
                self.U[:, i, j] = data[f"Total~Velocity~BLD_{i+1}~PAN_{j}~[m/s]"].to_numpy()
                self.alpha[:, i, j] = data[f"Angle~of~Attack~at~0.25c~BLD_{i+1}~PAN_{j}~[deg]"].to_numpy()

        # Convert the results from [deg] to [rad]
        self.alpha = np.radians(self.alpha)
        self.theta = np.radians(self.theta)
        self.phi = np.radians(self.phi)
        self.psi = np.radians(self.psi)
