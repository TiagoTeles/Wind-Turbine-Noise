"""
Author:   T. Moreira da Fonte Fonseca Teles
Email:    tmoreiradafont@tudelft.nl
Date:     2025-03-21
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
        displacement_thickness -- determine the displacement thickness

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

        for airfoil_name in os.listdir(airfoil_dir):
            self.airfoils.append(Airfoil(os.path.join(airfoil_dir, airfoil_name)))

        # Add the polar objects
        polar_dir = os.path.join(os.path.dirname(path), "Polars")

        for polar_name in os.listdir(polar_dir):
            self.polars.append(Polar(os.path.join(polar_dir, polar_name)))

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
        Determine the airfoil shape at a given position.

        Arguments:
            pos : float -- position along the blade, [m]

        Returns:
            airfoil : Airfoil -- interpolated airfoil
        """

        # Determine the neighbouring airfoils
        index_1 = self.data[self.data["pos"] <= pos].index.max()
        index_2 = self.data[self.data["pos"] >= pos].index.min()

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
        airfoil_path_1 = os.path.join("Airfoils", os.path.basename(airfoil_path_1))
        airfoil_path_2 = os.path.join("Airfoils", os.path.basename(airfoil_path_2))

        # Determine the interpolation fraction
        pos_1 = self.data["pos"].iloc[index_1]
        pos_2 = self.data["pos"].iloc[index_2]

        fraction = np.interp(pos, [pos_1, pos_2], [0, 1])

        # Interpolate the airfoil
        xfoil = XFoil(XFOIL_PATH, os.path.dirname(self.path))
        path_out = xfoil.interpolate(airfoil_path_1, airfoil_path_2, fraction)

        return Airfoil(path_out)


    def displacement_thickness(self, pos, Re, M, alpha, cutoff, probe_top, probe_bot, it):
        """
        Determine the boundary layer displacement thickness.

        Arguments:
            pos : float -- position along the blade, [-]
            Re : float -- Reynolds number, [-]
            M : float -- Mach number, [-]
            alpha : float -- angle of attack, [rad]
            cutoff : float -- radial cutoff, [-]
            probe_top : np.array -- probe location at the top surface, [-]
            probe_bot : np.array -- probe location at the bottom surface, [-]
            it : int -- number of iterations for XFoil, [-]

        Returns:
            delta_star_top : float -- top boundary layer displacement thickness, [m]
            delta_star_bot : float-- bottom boundary layer displacement thickness, [m]
        """

        airfoil = self.interpolate(pos)

        # Use a cutoff to ensure that XFoil converges
        # For cutoff=0.4, SPL_cutoff = SPL_tip - 19.90 dB
        if pos / self.data["pos"].iloc[-1] > cutoff:

            # Determine the path to the airfoil
            cwd = os.path.dirname(self.path)
            path = os.path.relpath(airfoil.path, cwd)

            # Run XFoil
            xfoil = XFoil(XFOIL_PATH, cwd)
            top, bot = xfoil.run(path, Re, M, alpha, it=it)

            # Determine the displacement thickness at the probe locations
            delta_star_top = np.interp(probe_top, top["x/c"], top["delta_star"])
            delta_star_bot = np.interp(probe_bot, bot["x/c"], bot["delta_star"])

        else:

            # Set the displacement thickness to zero
            delta_star_top = 0
            delta_star_bot = 0

        return delta_star_top, delta_star_bot
