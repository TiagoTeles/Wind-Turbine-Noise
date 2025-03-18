"""
Author:   T. Moreira da Fonte Fonseca Teles
Email:    tmoreiradafont@tudelft.nl
Date:     2025-03-18
License:  GNU GPL 3.0

Store the data from .bld files.

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
from settings import XFOIL_PATH
from xfoil import XFoil


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
        interpolate -- determine the interpolated airfoil

    Attributes:
        airfoils : list -- list of airfoil objects
        attributes : dict -- dictionary of attributes
        data : pd.DataFrame -- pos [m], chord [m], twist [rad], offset_x [m], offset_y [m], p_axis [-], and polar_file
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

        self.airfoils = []
        self.attributes = {}
        self.path = path
        self.polars = []

        # Check if the file exists
        if not os.path.isfile(path):
            print(f"No file found at {path}!")
            sys.exit(1)

        # Read the file
        self.read()

        # Add the airfoil objects
        airfoil_dir = os.path.join(os.path.dirname(path), "Airfoils")

        for f in os.listdir(airfoil_dir):
            self.airfoils.append(Airfoil(os.path.join(airfoil_dir, f)))

        if not airfoil_dir:
            print(f"No Airfoils found in {airfoil_dir}!")
            sys.exit(1)

        # Add the polar objects
        polar_dir = os.path.join(os.path.dirname(path), "Polars")

        for f in os.listdir(polar_dir):
            self.polars.append(Polar(os.path.join(polar_dir, f)))

        if not polar_dir:
            print(f"No Polars found in {polar_dir}!")
            sys.exit(1)

        # Precompute the interpolated airfoils
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

        # Open the file
        f = open(self.path, "r", encoding="utf-8")

        # Read the attributes
        for key, value in BLADE_DICT.items():
            self.attributes[key] = read(f, key, value["type"])

        # Read the pos, chord, twist, offset_x, offset_y, p_axis, and polar_file
        f.seek(0)

        self.data = pd.read_csv(f, names=["pos", "chord", "twist", "offset_x", "offset_y", \
                                          "p_axis", "polar_file"], skiprows=16, delimiter=r"\s+")

        # Format the twist and polar_file
        self.data["twist"] = np.radians(self.data["twist"])
        self.data["polar_file"] = self.data["polar_file"].str.replace("/","\\")

        # Close the file
        f.close()

    def interpolate(self, pos):
        """
        Determine the interpolated airfoil.

        Arguments:
            pos : float -- position along the blade, [m]
        
        Returns:
            airfoil : Airfoil -- interpolated airfoil
        """

        # Create the output directory
        dir_out = os.path.join(os.path.dirname(self.path), "Interpolated")

        if not os.path.exists(dir_out):
            os.makedirs(dir_out)

        # Determine the neighbouring airfoils
        index_1 = self.data[self.data["pos"] <= pos].index.max()
        index_2 = self.data[self.data["pos"] >= pos].index.min()

        r_1 = self.data["pos"].iloc[index_1]
        r_2 = self.data["pos"].iloc[index_2]

        polar_path_1  = self.data["polar_file"].iloc[index_1]
        polar_path_2 = self.data["polar_file"].iloc[index_2]

        airfoil_path_1  = None
        airfoil_path_2 = None

        for polar in self.polars:

            if os.path.basename(polar_path_1) == os.path.basename(polar.path):
                airfoil_path_1 = polar.attributes["FOILNAME"]

            if os.path.basename(polar_path_2) ==  os.path.basename(polar.path):
                airfoil_path_2 = polar.attributes["FOILNAME"]

        if airfoil_path_1 is None or airfoil_path_2 is None:
            print("Polar not found!")
            sys.exit(1)

        # Determine the absolute paths
        cwd = os.path.dirname(self.path)

        airfoil_path_1 = os.path.join("Airfoils", os.path.basename(airfoil_path_1))
        airfoil_path_2 = os.path.join("Airfoils", os.path.basename(airfoil_path_2))
        airfoil_path_interpolated = os.path.join("Interpolated", f"{pos}.afl")
        
        # Determine the interpolation fraction
        fraction = np.interp(pos, [r_1, r_2], [0, 1])

        # Interpolate the airfoil
        xfoil = XFoil(XFOIL_PATH, cwd)
        xfoil.interpolate(airfoil_path_1, airfoil_path_2, airfoil_path_interpolated, fraction)
        
        return Airfoil(os.path.join(os.path.dirname(self.path), airfoil_path_interpolated))
