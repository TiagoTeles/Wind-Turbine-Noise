"""
Author:   T. Moreira da Fonte Fonseca Teles
Email:    tmoreiradafont@tudelft.nl
Date:     2025-09-19
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

from source.QBlade.airfoil import Airfoil


class Polar:
    """
    A class to store the polar data.

    Methods:
        __init__ -- initialise the Polar class
        read -- read the .plr file

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

        # Add the Airfoil object
        airfoil_path = os.path.normpath(os.path.join(os.path.dirname(self.path), lines[8].split()[0]))
        self.airfoil = Airfoil(airfoil_path)

        # Close the file
        f.close()
