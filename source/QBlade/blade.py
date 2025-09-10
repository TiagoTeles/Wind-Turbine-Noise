"""
Author:   T. Moreira da Fonte Fonseca Teles
Email:    tmoreiradafont@tudelft.nl
Date:     2025-09-10
License:  GNU GPL 3.0

Store the blade data.

Classes:
    Blade

Functions:
    None

Exceptions:
    None
"""

import os
import sys

import numpy as np
import pandas as pd

from source.QBlade.polar import Polar


class Blade():
    """
    A class to store the blade data.

    Methods:
        __init__ -- initialise the Blade class
        read -- read the .bld file

    Attributes:
        path : str -- path to the .bld file
        geometry : pd.DataFrame -- blade geometry
    """

    def __init__(self, path):
        """
        Initialise the Blade class.

        Parameters:
            path : str -- path to the .bld file

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
        Read the .bld file.

        Parameters:
            None

        Returns:
            None
        """

        # Read the blade geometry
        self.geometry = pd.read_csv(self.path, delimiter=r"\s+", names=["radius", "chord", "twist", \
                                    "offset_x", "offset_y", "pitch_axis", "polar_path"], skiprows=16)

        # Convert the twist from [deg] to [rad]
        self.geometry["twist"] = np.radians(self.geometry["twist"])

        # Assign a Polar to each station
        self.geometry["polar"] = None

        for index, row in self.geometry.iterrows():

            # Determine the polar path
            polar_path = os.path.join(os.path.dirname(self.path), row["polar_path"])

            # Initialise the Polar
            self.geometry.at[index, "polar"] = Polar(polar_path)
