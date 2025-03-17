"""
Author:   T. Moreira da Fonte Fonseca Teles
Email:    tmoreiradafont@tudelft.nl
Date:     2025-03-11
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

import numpy as np
import pandas as pd


class Airfoil:
    """
    A class to store the airfoil data. 

    Methods:
        __init__ -- initialise the airfoil class
        read -- read the .afl file
        thickness -- determine the thickness of the airfoil

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

        # Check if the file exists
        if os.path.isfile(path):
            self.path = path
        else:
            print(f"No file found at {path}!")
            sys.exit(1)

        # Read the file
        self.read()

    def read(self):
        """
        Read the .afl file.

        Arguments:
            None

        Returns:
            None
        """

        # Open the file
        f = open(self.path, "r", encoding="utf-8")

        # Read the attributes
        self.attributes["AIRFOILNAME"] = f.readline().strip("\n")

        # Read x/c and y/c
        self.data = pd.read_csv(f, names=["x/c", "y/c"], skiprows=1, delimiter=r"\s+")

        # Close file
        f.close()

    def thickness(self, x):
        """
        Determine the thickness of the airfoil.

        Arguments:
            x : float -- x/c, [-]

        Returns:
            tc : float -- thickness of the airfoil
        """

        # Seperate the airfoil into top and bottom surfaces
        le_index = self.data["x/c"].idxmin()
        top = self.data.iloc[:le_index + 1][::-1]
        bot = self.data.iloc[le_index:]

        # Calculate the thickness
        y_top = np.interp(x, top["x/c"], top["y/c"])
        y_bot = np.interp(x, bot["x/c"], bot["y/c"])

        return y_top - y_bot
