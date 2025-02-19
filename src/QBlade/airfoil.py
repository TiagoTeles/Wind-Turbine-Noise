"""
Author:   T. Moreira da Fonte Fonseca Teles
Email:    tmoreiradafont@tudelft.nl
Date:     2025-02-19
License:  GNU GPL 3.0

Store data from .afl files.

Classes:
    Airfoil 

Functions:
    None

Exceptions:
    None
"""

import os
import sys

import pandas as pd


class Airfoil:
    """
    A class to store the airfoil data. 

    Methods:
        __init__ -- initialise the airfoil class
        read -- read the .afl file
        write -- write to the .afl file

    Attributes:
        attributes : dictionary -- dictionary of attributes
        data : pd.DataFrame -- x/c [-], y/c [-]
        path : str -- path to the .afl file
    """

    def __init__(self, path):
        """
        Initialise the Airfoil class.

        Arguments:
            path : str -- path to the .afl file

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
        Read the .afl file.

        Arguments:
            None

        Returns:
            None
        """

        # Open file
        f = open(self.path, "r", encoding="utf-8")

        # Read attributes
        self.attributes["AIRFOILNAME"] = f.readline().strip("\n")

        # Read x/c and y/c
        f.seek(0)

        self.data = pd.read_csv(f, names=["x/c", "y/c"], skiprows=1, delimiter=r"\s+")

        # Close file
        f.close()

    def write(self):
        pass

# if __name__ == "__main__":

#     # Create an Airfoil instance
#     airfoil = Airfoil("data\\turbines\\DTU_10MW\\Aero\\Airfoils\\FFA_W3_241.afl")

#     # Print attributes
#     print(f"AIRFOILNAME: {airfoil.attributes['AIRFOILNAME']}")

#     # Print data
#     print(airfoil.data)
