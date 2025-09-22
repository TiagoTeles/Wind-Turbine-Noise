"""
Author:   T. Moreira da Fonte Fonseca Teles
Email:    tmoreiradafont@tudelft.nl
Date:     2025-09-22
License:  GNU GPL 3.0

Store the simulation data.

Classes:
    Simulation

Functions:
    None

Exceptions:
    None
"""

import os

from source.QBlade.turbine import Turbine


class Simulation:
    """
    A class to store the simulation data.

    Methods:
        __init__ -- initialise the Simulation class

    Attributes:
        path : str -- path to the .sim file
        turbine: Turbine -- turbine object
    """

    def __init__(self, path):
        """
        Initialise the Simulation class.

        Parameters:
            path : str -- path to the .sim file

        Returns:
            None
        """

        self.path = path

        # Open the file
        f = open(self.path, "r", encoding="utf-8")
        lines = f.readlines()

        # Add the Turbine object
        turbine_path = os.path.normpath(os.path.join(os.path.dirname(self.path), lines[16].split()[0]))
        self.turbine = Turbine(turbine_path)

        # Close the file
        f.close()
