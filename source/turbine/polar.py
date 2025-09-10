"""
Author:   T. Moreira da Fonte Fonseca Teles
Email:    tmoreiradafont@tudelft.nl
Date:     2025-09-10
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

from source.turbine.airfoil import Airfoil


class Polar:
    """
    A class to store the polar data.

    Methods:
        __init__ -- initialise the polar class
        read -- read the .plr file
    
    Attributes:
        path : str -- path to the .plr file
        airfoil : Airfoil -- airfoil object
        reynolds : float -- Reynolds number
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

        # Read the file
        lines = f.readlines()

        # Determine the airfoil path
        airfoil_path = lines[8].split()[0]
        polar_path = os.path.dirname(self.path)
        path = os.path.normpath(os.path.join(polar_path, airfoil_path))

        # Initialise the airfoil
        self.airfoil = Airfoil(path)

        # Read the Reynolds number
        self.reynolds = float(lines[13].split()[1])

        # Close the file
        f.close()

        # Read alpha, Cl, Cd, and Cm
        self.coefficients = pd.read_csv(self.path, delimiter=r"\s+", \
                                        names=["alpha", "Cl", "Cd", "Cm"], skiprows=17)

        # Convert alpha from [deg] to [rad]
        self.coefficients["alpha"] = np.radians(self.coefficients["alpha"])
