"""
Author:   T. Moreira da Fonte Fonseca Teles
Email:    tmoreiradafont@tudelft.nl
Date:     2025-05-15
License:  GNU GPL 3.0

Store the data from .plr files.

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

from QBlade.io import read


POLAR_DICT = {
    "POLARNAME":    {"type":   str}, # Polar name
    "FOILNAME":     {"type":   str}, # Path to the .afl file
    "THICKNESS":    {"type": float}, # Airfoil thickness, [-]
    "ISDECOMPOSED": {"type":  bool}, # Is the polar is decomposed?
    }

class Polar:
    """
    A class to store the polar data.

    Methods:
        __init__ -- initialise the polar class
        read -- read the .plr file

    Attributes:
        attributes : dict -- dictionary of attributes
        data : pd.DataFrame -- alpha [rad], CL [-], CD [-], and CM [-]
        path : str -- path to the .plr file
    """

    def __init__(self, path):
        """
        Initialise the Polar class.

        Arguments:
            path : str -- path to the .plr file

        Returns:
            None
        """

        self.attributes = {}
        self.path = path

        # Check if the file exists
        if not os.path.isfile(path):
            print(f"No file found at {path}!")
            sys.exit(1)

        # Read file
        self.read()

    def read(self):
        """
        Read the .plr file.

        Arguments:
            None

        Returns:
            None
        """

        # Open the file
        f = open(self.path, "r", encoding="utf-8")

        # Read the attributes
        for key, value in POLAR_DICT.items():
            self.attributes[key] = read(f, key, value["type"])

        # Format the attributes
        self.attributes["FOILNAME"] = self.attributes["FOILNAME"].replace("/", "\\")
        self.attributes["THICKNESS"] /= 100

        # Read the AoA, CL, CD, and CM
        f.seek(0)
        self.data = pd.read_csv(f, names=["alpha", "Cl", "Cd", "Cm"], skiprows=17, delimiter=r"\s+")

        # Format the AoA
        self.data["alpha"] = np.radians(self.data["alpha"])

        # Close file
        f.close()
