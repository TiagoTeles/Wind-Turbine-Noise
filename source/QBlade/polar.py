"""
Author:   T. Moreira da Fonte Fonseca Teles
Email:    tmoreiradafont@tudelft.nl
Date:     2025-09-17
License:  GNU GPL 3.0

Store the polar data.

Classes:
    Polar

Functions:
    None

Exceptions:
    None
"""

import os
import sys

import numpy as np
import pandas as pd

from source.QBlade.airfoil import Airfoil


class Polar:
    """
    A class to store the polar data.

    Methods:
        __init__ -- initialise the Polar class
        read -- read the .plr file

    Attributes:
        path : str -- path to the .plr file
        reynolds : float -- Reynolds number
        airfoil : Airfoil -- airfoil object
        coefficients : pd.DataFrame -- alpha, Cl, Cd, and Cm
    """

    def __init__(self, path):
        """
        Initialise the Polar class.

        Parameters:
            path : str -- path to the .plr file

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
        Read the .plr file.

        Parameters:
            None

        Returns:
            None
        """

        # Open the file
        f = open(self.path, "r", encoding="utf-8")
        lines = f.readlines()

        # Read the Reynolds number
        self.reynolds = float(lines[13].split()[1])

        # Add the Airfoil object
        airfoil_path = os.path.normpath(os.path.join(os.path.dirname(self.path), lines[8].split()[0]))
        self.airfoil = Airfoil(airfoil_path)

        # Close the file
        f.close()

        # Read the polar coefficients
        self.coefficients = pd.read_csv(self.path, delimiter=r"\s+", \
                                        names=["alpha", "Cl", "Cd", "Cm"], skiprows=17)

        # Convert alpha from [deg] to [rad]
        self.coefficients["alpha"] = np.radians(self.coefficients["alpha"])
