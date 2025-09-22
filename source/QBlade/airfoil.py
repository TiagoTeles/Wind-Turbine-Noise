"""
Author:   T. Moreira da Fonte Fonseca Teles
Email:    tmoreiradafont@tudelft.nl
Date:     2025-09-22
License:  GNU GPL 3.0

Store the airfoil data.

Classes:
    Airfoil

Functions:
    None

Exceptions:
    None
"""

import numpy as np
import pandas as pd


class Airfoil:
    """
    A class to store the airfoil data.

    Methods:
        __init__ -- initialise the Airfoil class
        thickness -- determine the airfoil thickness

    Attributes:
        path : str -- path to the .afl file
        coordinates : pd.DataFrame -- airfoil coordinates
    """

    def __init__(self, path):
        """
        Initialise the Airfoil class.

        Parameters:
            path : str -- path to the .afl file

        Returns:
            None
        """

        self.path = path

        # Read the airfoil coordinates
        self.coordinates = pd.read_csv(self.path, delimiter=r"\s+", names=["x/c", "y/c"], skiprows=1)

    def thickness(self, x_c):
        """
        Determine the airfoil thickness.

        Parameters:
            x_c : float -- x/c, [-]

        Returns:
            t_c : float -- t/c, [-]
        """

        # Separate the upper and lower airfoil surfaces
        index = self.coordinates["x/c"].idxmin()
        upper = self.coordinates.iloc[:index + 1][::-1]
        lower = self.coordinates.iloc[index:]

        # Determine the airfoil thickness
        y_c_upper = np.interp(x_c, upper["x/c"], upper["y/c"])
        y_c_lower = np.interp(x_c, lower["x/c"], lower["y/c"])

        return y_c_upper - y_c_lower
