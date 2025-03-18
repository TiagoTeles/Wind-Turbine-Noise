"""
Author:   T. Moreira da Fonte Fonseca Teles
Email:    tmoreiradafont@tudelft.nl
Date:     2025-03-18
License:  GNU GPL 3.0

Store the data from .trb files.

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

from QBlade.blade import Blade
from QBlade.misc import read, write


TURBINE_DICT = {
    "OBJECTNAME":       {"type":   str}, # Name of the turbine object
    "BLADEFILE":        {"type":   str}, # Path of the blade file
    "TURBTYPE":         {"type":   int}, # Turbine type (0=HAWT or 1=VAWT)
    "NUMBLADES":        {"type":   int}, # Number of blades, [-]
    "ROTORCONFIG":      {"type":   int}, # Rotor configuration (0=UPWIND or 1=DOWNWIND)
    "ROTATIONALDIR":    {"type":   int}, # Direction of rotor rotation (0=STANDARD or 1=REVERSED)
    "DISCTYPE":         {"type":   int}, # Type of rotor discretization (0=from Bladetable, 1=linear, 2=cosine)
    "INTPTYPE":         {"type":   int}, # Type of rotor interpolation (0=linear, 1=Akima splines)
    "NUMPANELS":        {"type":   int}, # Number of aerodynamic panels per blade, [-]
    "OVERHANG":         {"type": float}, # Rotor overhang, [m]
    "SHAFTTILT":        {"type": float}, # Shaft tilt angle, [rad]
    "ROTORCONE":        {"type": float}, # Rotor cone angle, [rad]
    "CLEARANCE":        {"type": float}, # Rotor clearance to ground, [m]
    "XTILT":            {"type": float}, # Rotor x-tilt angle, [rad]
    "YTILT":            {"type": float}, # Rotor y-tilt angle, [rad]
    "TOWERHEIGHT":      {"type": float}, # Tower height, [m]
    "TOWERTOPRAD":      {"type": float}, # Tower top radius, [m]
    "TOWERBOTRAD":      {"type": float}, # Tower bottom radius, [m]
    "DYNSTALLTYPE":     {"type":   int}, # Dynamic stall model (0=none, 1=OYE, 2=IAG, 3=GORMONT-BERG, 4=ATEFLAP)
    "TF_OYE":           {"type": float}, # Tf constant for the OYE dynamic stall model, [-]
    "AM_GB":            {"type": float}, # Am constant for the GORMONT-BERG dynamic stall model, [-]
    "TF_ATE":           {"type": float}, # Tf constant for the ATEFLAP dynamic stall model, [-]
    "TP_ATE":           {"type": float}, # Tp constant for the ATEFLAP dynamic stall model, [-]
    "UNSTEADYAERO":     {"type":  bool}, # Include unsteady non-circulatory aerodynamics?
    "2PLIFTDRAG":       {"type":  bool}, # Include the 2 point lift drag correction?
    "HIMMELSKAMP":      {"type":  bool}, # Include the Himmelskamp Stall delay? 
    "TOWERSHADOW":      {"type":  bool}, # Include the tower shadow effect?
    "TOWERDRAG":        {"type": float}, # Tower drag coefficient, [-]
    "WAKETYPE":         {"type":   int}, # Wake type (0=free vortex wake, 1=unsteady BEM)
    "WAKEINTTYPE":      {"type":   int}, # Wake integration type (0=EF, 1=ET, 2=PC, 3=PC2B)
    "WAKEROLLUP":       {"type":  bool}, # Calculate wake self-induction?
    "TRAILINGVORT":     {"type":  bool}, # Include trailing vortex elements?
    "SHEDVORT":         {"type":  bool}, # Include shed vortex elements?
    "CONVECTIONTYPE":   {"type":   int}, # Wake convection type (0=BL, 1=HH, 2=LOC)
    "WAKERELAXATION":   {"type": float}, # Wake relaxation factor, [-]
    "FIRSTWAKEROW":     {"type": float}, # First wake row length, [-]
    "MAXWAKESIZE":      {"type":   int}, # Maximum number of wake elements, [-]
    "MAXWAKEDIST":      {"type": float}, # Maximum wake distance from the rotor plane, [-]
    "WAKEREDUCTION":    {"type": float}, # Wake reduction factor, [-]
    "WAKELENGTHTYPE":   {"type":   int}, # Wake length type (0=counted in rotor revolutions, 1=counted in time steps)
    "CONVERSIONLENGTH": {"type": float}, # Wake conversion length, [-]
    "NEARWAKELENGTH":   {"type": float}, # Near wake length, [-]
    "ZONE1LENGTH":      {"type": float}, # Wake zone 1 length, [-]
    "ZONE2LENGTH":      {"type": float}, # Wake zone 2 length, [-]
    "ZONE3LENGTH":      {"type": float}, # Wake zone 3 length, [-]
    "ZONE1FACTOR":      {"type":   int}, # Wake zone 1 streamwise factor, [-]
    "ZONE2FACTOR":      {"type":   int}, # Wake zone 2 streamwise factor, [-]
    "ZONE3FACTOR":      {"type":   int}, # Wake zone 3 streamwise factor, [-]
    "ZONE1FACTOR_S":    {"type":   int}, # Wake zone 1 spanwise factor, [-]
    "ZONE2FACTOR_S":    {"type":   int}, # Wake zone 2 spanwise factor, [-]
    "ZONE3FACTOR_S":    {"type":   int}, # Wake zone 3 spanwise factor, [-]
    "BOUNDCORERADIUS":  {"type": float}, # Fixed core radius of the bound blade vortex, [-]
    "WAKECORERADIUS":   {"type": float}, # Initial core radius of the free wake vortex, [-]
    "VORTEXVISCOSITY":  {"type": float}, # Turbulent vortex viscosity, [-]
    "VORTEXSTRAIN":     {"type":  bool}, # Calculate vortex strain?
    "MAXSTRAIN":        {"type":   int}, # Maximum element strain, [-]
    "GAMMARELAXATION":  {"type": float}, # Relaxation factor used in the gamma iteration, [-]
    "GAMMAEPSILON":     {"type": float}, # Relative gamma convergence criteria, [-]
    "GAMMAITERATIONS":  {"type":   int}, # Maximum number of gamma iterations, [-]
    "POLARDISC":        {"type":   int}, # Polar discretization for the unsteady BEM, [-]
    "BEMTIPLOSS":       {"type":  bool}, # Use BEM tip loss factor?
    "BEMSPEEDUP":       {"type": float}, # Initial BEM convergence acceleration time, [s]
    "STRUCTURALFILE":   {"type":   str}, # Input file for the structural model (leave blank if unused)
    "GEOMSTIFFNESS":    {"type":  bool}, # Enable geometric stiffness?
    "AEROPANELLOADS":   {"type":  bool}, # Enable distributed aero panel loads and gradients?
    "CONTROLLERTYPE":   {"type":   int}, # Type of turbine controller (0=none, 1=BLADED, 2=DTU, 3=TUB)
    "CONTROLLERFILE":   {"type":   str}, # Controller file name (leave blank if unused)
    "PARAMETERFILE":    {"type":   str}, # Controller parameter file name (leave blank if unused)
    }

class Turbine:
    """
    A class to store the turbine data.

    Methods:
        __init__ -- initialise the turbine class
        read -- read the .trb file
        write -- write to the .trb file

    Attributes:
        attributes : dict -- dictionary of attributes
        blade: Blade -- blade object
        path : str -- path to the .trb file
    """

    def __init__(self, path):
        """
        Initialise the Turbine class.

        Arguments:
            path : str -- path to the .trb file

        Returns:
            None
        """

        self.attributes = {}
        self.path = path

        # Check if the file exists
        if not os.path.isfile(path):
            print(f"No file found at {path}!")
            sys.exit(1)

        # Read the file
        self.read()

        # Add the blade object
        blade_path = os.path.join(os.path.dirname(self.path), self.attributes["BLADEFILE"])
        self.blade = Blade(blade_path)

    def read(self):
        """
        Read the .trb file.

        Arguments:
            None

        Returns:
            None
        """

        # Open the file
        f = open(self.path, "r", encoding="utf-8")

        # Read the attributes
        for key, value in TURBINE_DICT.items():
            self.attributes[key] = read(f, key, value["type"])

        # Format the attributes
        self.attributes["BLADEFILE"]      = self.attributes["BLADEFILE"].replace("/", "\\")
        self.attributes["STRUCTURALFILE"] = self.attributes["STRUCTURALFILE"].replace("/", "\\")
        self.attributes["CONTROLLERFILE"] = self.attributes["CONTROLLERFILE"].replace("/", "\\")
        self.attributes["PARAMETERFILE"]  = self.attributes["PARAMETERFILE"].replace("/", "\\")

        self.attributes["SHAFTTILT"] = np.radians(self.attributes["SHAFTTILT"])
        self.attributes["ROTORCONE"] = np.radians(self.attributes["ROTORCONE"])
        self.attributes["XTILT"]     = np.radians(self.attributes["XTILT"])
        self.attributes["YTILT"]     = np.radians(self.attributes["YTILT"])

        # Close the file
        f.close()

    def write(self, key, value):

        # Set the value in the attributes dictionary
        self.attributes[key] = value

        # Set the value in .trb file
        f = open(self.path, "r+", encoding="utf-8")
        write(f, key, value)
        f.close()
