"""
Author:   T. Moreira da Fonte Fonseca Teles
Email:    tmoreiradafont@tudelft.nl
Date:     2025-02-19
License:  GNU GPL 3.0

Store data from .plr files.

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

from misc import read


POLAR_DICT = {
    "POLARNAME":    {"type":   str}, # Polar name
    "FOILNAME":     {"type":   str}, # Path to the .afl file
    "THICKNESS":    {"type": float}, # Airfoil thickness, [-]
    "ISDECOMPOSED": {"type":  bool}, # Is the polar is decomposed?
    "REYNOLDS":     {"type": float}, # Reynolds number, [-]
    }

class Polar:
    """
    A class to store the polar data.

    Methods:
        __init__ -- initialise the polar class
        read -- read the .plr file
        write -- write to the .plr file

    Attributes:
        attributes : dict -- dictionary of attributes
        data : pd.DataFrame -- AOA [rad], CL [-], CD [-], and CM [-]
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

        # Check if file exists
        if os.path.isfile(path):
            self.path = path
        else:
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

        # Open file
        f = open(self.path, "r", encoding="utf-8")

        # Read attributes
        for key, value in POLAR_DICT.items():
            self.attributes[key] = read(f, key, value["type"])

        self.attributes["FOILNAME"] = self.attributes["FOILNAME"].replace("/", "\\")
        self.attributes["THICKNESS"] /= 100

        # Read AOA, CL, CD, and CM
        f.seek(0)

        self.data = pd.read_csv(f, names=["AoA", "Cl", "Cd", "Cm"], skiprows=17, delimiter=r"\s+")
        self.data["AoA"] = np.radians(self.data["AoA"])

        # Close file
        f.close()

    def write(self):
        pass

# if __name__ == "__main__":

#     # Create a Polar instance
#     polar = Polar("data\\turbines\\DTU_10MW\\Aero\\Polars\\FFA_W3_241_t24.1_dtu_10mw_Polar_RE1.00E+06.plr")

#     # Print attributes
#     for key, value in polar.attributes.items():
#         print(f"{key}: {value}")

#     # Print data
#     print(polar.data)
