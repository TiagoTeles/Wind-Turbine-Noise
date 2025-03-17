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

from QBlade.airfoil import Airfoil
from QBlade.misc import read
from QBlade.polar import Polar
from xfoil import XFoil


BLADE_DICT = {
    "OBJECTNAME":    {"type":  str}, # Name of the blade object
    "ROTORTYPE":     {"type":  str}, # Rotor type
    "INVERTEDFOILS": {"type": bool}, # Invert the airfoils?
    "NUMBLADES":     {"type":  int}, # Number of blades, [-]
    }

XFOIL_PATH = "bin\\XFoil\\xfoil.exe"

class Blade():
    """
    A class to store the aerodynamic data.

    Methods:
        __init__ -- initialise the blade class
        read -- read the .bld file
        write -- write to the .bld file

    Attributes:
        airfoils : list -- list of airfoil objects
        attributes : dict -- dictionary of attributes
        data : pd.DataFrame -- pos [m], chord [m], twist [rad], offset_x [m], offset_y [m], p_axis [-], and polar_path
        path : str -- path to the .bld file
        polars : list -- list of polar objects
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

        # Precompute interpolated airfoils
        self.data.insert(6, "airfoil", None)

        for i, pos in enumerate(self.data["pos"]):
            self.data.loc[i, "airfoil"] = self.interpolate(pos)

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

        self.data["polar_path"] = self.data["polar_path"].str.replace("/","\\")
        self.data["twist"] = np.radians(self.data["twist"])

        # Close file
        f.close()

    def write(self):
        pass

    def interpolate(self, pos):

        if not os.path.exists(os.path.join(os.path.dirname(self.path), "Interpolated")):
            os.makedirs(os.path.join(os.path.dirname(self.path), "Interpolated"))

        # Determine the neighbouring airfoils
        index_inboard  = self.data[self.data["pos"] <= pos].index.max()
        index_outboard = self.data[self.data["pos"] >= pos].index.min()

        r_inboard  = self.data["pos"].iloc[index_inboard]
        r_outboard = self.data["pos"].iloc[index_outboard]

        polar_path_inboard  = self.data["polar_path"].iloc[index_inboard]
        polar_path_outboard = self.data["polar_path"].iloc[index_outboard]

        airfoil_path_inboard  = None
        airfoil_path_outboard = None

        for polar in self.polars:

            if os.path.basename(polar_path_inboard) == os.path.basename(polar.path):
                airfoil_path_inboard = polar.attributes["FOILNAME"]

            if os.path.basename(polar_path_outboard) ==  os.path.basename(polar.path):
                airfoil_path_outboard = polar.attributes["FOILNAME"]

        if airfoil_path_inboard is None or airfoil_path_outboard is None:
            print("Polar not found!")
            sys.exit(1)

        airfoil_path_inboard = os.path.relpath(os.path.join(os.path.dirname(polar_path_inboard), airfoil_path_inboard))
        airfoil_path_outboard =  os.path.relpath(os.path.join(os.path.dirname(polar_path_outboard), airfoil_path_outboard))
        airfoil_path_interpolated = os.path.join("Interpolated", f"{pos}.afl")
        
        # Interpolate the airfoil
        xfoil = XFoil(XFOIL_PATH, os.path.dirname(self.path))
        xfoil.interpolate(airfoil_path_inboard, airfoil_path_outboard, airfoil_path_interpolated, np.interp(pos, [r_inboard, r_outboard], [0, 1]))
        
        return Airfoil(os.path.join(os.path.dirname(self.path), airfoil_path_interpolated))