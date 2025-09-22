"""
Author:   T. Moreira da Fonte Fonseca Teles
Email:    tmoreiradafont@tudelft.nl
Date:     2025-09-22
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

from source.QBlade.airfoil import Airfoil


class Polar:
    """
    A class to store the polar data.

    Methods:
        __init__ -- initialise the Polar class

    Attributes:
        path : str -- path to the .plr file
        airfoil : Airfoil -- airfoil object
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

        # Read the file
        f = open(self.path, "r", encoding="utf-8")
        lines = f.readlines()

        # Add the Airfoil object
        airfoil_path = os.path.normpath(os.path.join(os.path.dirname(self.path), lines[8].split()[0]))
        self.airfoil = Airfoil(airfoil_path)

        # Close the file
        f.close()
