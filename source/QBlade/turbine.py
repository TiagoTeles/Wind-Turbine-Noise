"""
Author:   T. Moreira da Fonte Fonseca Teles
Email:    tmoreiradafont@tudelft.nl
Date:     2025-09-10
License:  GNU GPL 3.0

Store the turbine data.

Classes:
    Turbine

Functions:
    None

Exceptions:
    None
"""

import os
import sys

import numpy as np

from source.QBlade.blade import Blade


class Turbine:
    """
    A class to store the turbine data.

    Methods:
        __init__ -- initialise the Turbine class
        read -- read the .trb file

    Attributes:
        path : str -- path to the .trb file
        n_blades : int -- number of blades
        n_panels : int -- number of panels
        overhang : float -- overhang, [m]
        shaft_tilt : float -- shaft tilt angle, [rad]
        rotor_cone : float -- rotor cone angle, [rad]
        tower_height : float -- tower height, [m]
        blade: Blade -- blade object
    """

    def __init__(self, path):
        """
        Initialise the Turbine class.

        Parameters:
            path : str -- path to the .trb file

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
        Read the .trb file.

        Parameters:
            None

        Returns:
            None
        """

        # Open the file
        f = open(self.path, "r", encoding="utf-8")
        lines = f.readlines()

        # Read the turbine geometry
        self.n_blades = int(lines[12].split()[0])
        self.n_panels = int(lines[17].split()[0])
        self.overhang = float(lines[21].split()[0])
        self.shaft_tilt = float(lines[22].split()[0])
        self.rotor_cone = float(lines[23].split()[0])
        self.tower_height = float(lines[27].split()[0])

        # Convert the angles cone from [deg] to [rad]
        self.shaft_tilt = np.radians(self.shaft_tilt)
        self.rotor_cone = np.radians(self.rotor_cone)

        # Add the blade object
        blade_path = os.path.join(os.path.dirname(self.path), lines[10].split()[0])
        self.blade = Blade(blade_path)

        # Close the file
        f.close()
