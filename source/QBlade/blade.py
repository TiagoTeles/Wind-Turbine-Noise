"""
Author:   T. Moreira da Fonte Fonseca Teles
Email:    tmoreiradafont@tudelft.nl
Date:     2025-05-15
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
from QBlade.io import read
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
        interpolate -- determine the airfoil at a given radius
        displacement_thickness -- determine the displacement thickness

    Attributes:
        airfoils : list -- list of airfoil objects
        attributes : dict -- dictionary of attributes
        data : pd.DataFrame -- radius [m], chord [m], twist [rad], offset_x [m], 
                               offset_y [m], pitch_axis [-], and polar_file
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

        # Read the radius, chord, twist, offset_x, offset_y, pitch_axis, and polar_file
        f.seek(0)

        self.data = pd.read_csv(f, names=["radius", "chord", "twist", "offset_x", "offset_y", \
                                          "pitch_axis", "polar_file"], skiprows=16, delimiter=r"\s+")

        # Format the twist and polar_file
        self.data["twist"] = np.radians(self.data["twist"])
        self.data["polar_file"] = self.data["polar_file"].str.replace("/","\\")

        # Close the file
        f.close()

    def interpolate(self, radius):
        """
        Determine the airfoil shape at a given radius.

        Arguments:
            radius : float -- radius along the blade, [m]

        Returns:
            airfoil : Airfoil -- interpolated airfoil
        """

        # Determine the neighbouring airfoils
        index_1 = self.data[self.data["radius"] <= radius].index.max()
        index_2 = self.data[self.data["radius"] >= radius].index.min()

        polar_path_1 = self.data["polar_file"].iloc[index_1]
        polar_path_2 = self.data["polar_file"].iloc[index_2]

        airfoil_path_1 = None
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
        radius_1 = self.data["radius"].iloc[index_1]
        radius_2 = self.data["radius"].iloc[index_2]

        fraction = np.interp(radius, [radius_1, radius_2], [0, 1])

        # Interpolate the airfoil
        xfoil = XFoil(XFOIL_PATH, os.path.dirname(self.path))
        path_out = xfoil.interpolate(airfoil_path_1, airfoil_path_2, fraction)

        return Airfoil(path_out)

    def displacement_thickness(self, radius, Re, M, alpha, cutoff, probe_top, probe_bot, x_tr_top, x_tr_bot, n_crit, max_iter):
        """
        Determine the boundary layer displacement thickness.

        Arguments:
            radius : float -- radius along the blade, [-]
            Re : float -- Reynolds number, [-]
            M : float -- Mach number, [-]
            alpha : float -- angle of attack, [rad]
            cutoff : float -- radial cutoff, [-]
            probe_top : np.array -- top probe position, [-]
            probe_bot : np.array -- bottom probe position, [-]
            x_tr_top : float -- top transition location, [-]
            x_tr_bot : float -- bottom transition location, [-]
            n_crit : float -- critical amplification factor, [-]
            max_iter : int -- maximum number of XFOIL iterations, [-]

        Returns:
            delta_star_top : float -- top boundary layer displacement thickness, [m]
            delta_star_bot : float-- bottom boundary layer displacement thickness, [m]
        """

        # Determine the airfoil
        airfoil = self.interpolate(radius)

        # Determine the chord
        chord = np.interp(radius, self.data["radius"], self.data["chord"])

        # Use a cutoff to ensure that XFoil converges
        if radius / self.data["radius"].iloc[-1] >= cutoff:

            # Determine the path to the airfoil
            cwd = os.path.dirname(self.path)
            path = os.path.relpath(airfoil.path, cwd)

            # Run XFoil
            xfoil = XFoil(XFOIL_PATH, cwd)
            top, bot = xfoil.run(path, Re, M, alpha, x_tr_bot, x_tr_bot, n_crit, max_iter)

            # Determine the displacement thickness at the probe locations
            delta_star_top = np.interp(probe_top, top["x/c"], top["delta_star/c"]) * chord
            delta_star_bot = np.interp(probe_bot, bot["x/c"], bot["delta_star/c"]) * chord

        else:

            # Set the displacement thickness to zero
            delta_star_top = 0
            delta_star_bot = 0

        return delta_star_top, delta_star_bot
