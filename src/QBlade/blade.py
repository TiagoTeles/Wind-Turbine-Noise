"""
Author:   T. Moreira da Fonte Fonseca Teles
Email:    tmoreiradafont@tudelft.nl
Date:     2025-02-19
License:  GNU GPL 3.0

Store data from .bld files.

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

from airfoil import Airfoil
from misc import read
from polar import Polar


BLADE_DICT = {
    "OBJECTNAME":    {"type":  str}, # Name of the blade object
    "ROTORTYPE":     {"type":  str}, # Rotor type
    "INVERTEDFOILS": {"type": bool}, # Invert the airfoils?
    "NUMBLADES":     {"type":  int}, # Number of blades, [-]
    }

class Blade():
    """
    A class to store the aerodynamic data.

    Methods:
        __init__ -- initialise the blade class
        read -- read the .bld file
        write -- write to the .bld file

    Attributes:
        airfoils : list -- list of Airfoil objects
        attributes : dict -- dictionary of attributes
        data : pd.DataFrame -- pos [m], chord [m], twist [rad], offset_x [m], offset_y [m], p_axis [-], and polar_path
        path : str -- path to the .bld file
        polars : list -- list of Polar objects
    """

    def __init__(self, path):
        """
        Initialise the Blade class.

        Arguments:
            path : str -- path to the .bld file

        Returns:
            None
        """

        self.attributes = {}

        self.airfoils = []
        self.polars = []

        # Check if file exists
        if os.path.isfile(path):
            self.path = path
        else:
            print(f"No file found at {path}!")
            sys.exit(1)

        # Read file
        self.read()

        # Add airfoil objects
        airfoil_dir = os.path.join(os.path.dirname(path), "Airfoils")

        if os.path.isdir(airfoil_dir):
            for f in os.listdir(airfoil_dir):
                self.airfoils.append(Airfoil(os.path.join(airfoil_dir, f)))
        else:
            print(f"No 'Airfoils' folder found at {airfoil_dir}!")
            sys.exit(1)

        # Add polar objects
        polar_dir = os.path.join(os.path.dirname(path), "Polars")

        if os.path.isdir(polar_dir):
            for f in os.listdir(polar_dir):
                self.polars.append(Polar(os.path.join(polar_dir, f)))
        else:
            print(f"No 'Polars' folder found at {polar_dir}!")
            sys.exit(1)

    def read(self):
        """
        Read the .bld file.

        Arguments:
            None

        Returns:
            None
        """

        # Open file
        f = open(self.path, "r", encoding="utf-8")

        # Read attributes
        for key, value in BLADE_DICT.items():
            self.attributes[key] = read(f, key, value["type"])

        # Read pos, chord, twist, offset_x, offset_y, p_axis, and polar_path
        f.seek(0)

        self.data = pd.read_csv(f, names=["pos", "chord", "twist", "offset_x", "offset_y", \
                                          "p_axis", "polar_path"], skiprows=16, delimiter=r"\s+")

        self.data["twist"] = np.radians(self.data["twist"])
        self.data["polar_path"] = os.path.dirname(self.path) + os.sep + self.data["polar_path"].str.replace("/","\\")

        # Close file
        f.close()

    def write(self):
        pass

# if __name__ == "__main__":

#     # Create a Blade instance
#     blade = Blade("data\\turbines\\DTU_10MW\\Aero\\DTU_10MW.bld")

#     # Print attributes
#     for key, value in blade.attributes.items():
#         print(f"{key}: {value}")

#     # Print data
#     print(blade.data)
