"""
Author:   T. Moreira da Fonte Fonseca Teles
Email:    tmoreiradafont@tudelft.nl
Date:     2025-02-19
License:  GNU GPL 3.0

Store data from .trb files.

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

from blade import Blade
from misc import read


TURBINE_DICT = {
    "OBJECTNAME":       {"type":   str, "pos": 0}, # Name of the turbine object
    "BLADEFILE":        {"type":   str, "pos": 0}, # Path of the blade file
    "TURBTYPE":         {"type":   int, "pos": 0}, # Turbine type (0=HAWT or 1=VAWT)
    "NUMBLADES":        {"type":   int, "pos": 0}, # Number of blades, [-]
    "ROTORCONFIG":      {"type":   int, "pos": 0}, # Rotor configuration (0=UPWIND or 1=DOWNWIND)
    "ROTATIONALDIR":    {"type":   int, "pos": 0}, # Direction of rotor rotation (0=STANDARD or 1=REVERSED)
    "DISCTYPE":         {"type":   int, "pos": 0}, # Type of rotor discretization (0=from Bladetable, 1=linear, 2=cosine)
    "INTPTYPE":         {"type":   int, "pos": 0}, # Type of rotor interpolation (0=linear, 1=Akima splines)
    "NUMPANELS":        {"type":   int, "pos": 0}, # Number of aerodynamic panels per blade, [-]
    "OVERHANG":         {"type": float, "pos": 0}, # Rotor overhang, [m]
    "SHAFTTILT":        {"type": float, "pos": 0}, # Shaft tilt angle, [rad]
    "ROTORCONE":        {"type": float, "pos": 0}, # Rotor cone angle, [rad]
    "CLEARANCE":        {"type": float, "pos": 0}, # Rotor clearance to ground, [m]
    "XTILT":            {"type": float, "pos": 0}, # Rotor x-tilt angle, [rad]
    "YTILT":            {"type": float, "pos": 0}, # Rotor y-tilt angle, [rad]
    "TOWERHEIGHT":      {"type": float, "pos": 0}, # Tower height, [m]
    "TOWERTOPRAD":      {"type": float, "pos": 0}, # Tower top radius, [m]
    "TOWERBOTRAD":      {"type": float, "pos": 0}, # Tower bottom radius, [m]
    "DYNSTALLTYPE":     {"type":   int, "pos": 0}, # Dynamic stall model (0=none, 1=OYE, 2=IAG, 3=GORMONT-BERG, 4=ATEFLAP)
    "TF_OYE":           {"type": float, "pos": 0}, # Tf constant for the OYE dynamic stall model, [-]
    "AM_GB":            {"type": float, "pos": 0}, # Am constant for the GORMONT-BERG dynamic stall model, [-]
    "TF_ATE":           {"type": float, "pos": 0}, # Tf constant for the ATEFLAP dynamic stall model, [-]
    "TP_ATE":           {"type": float, "pos": 0}, # Tp constant for the ATEFLAP dynamic stall model, [-]
    "UNSTEADYAERO":     {"type":  bool, "pos": 0}, # Include unsteady non-circulatory aerodynamics?
    "2PLIFTDRAG":       {"type":  bool, "pos": 0}, # Include the 2 point lift drag correction?
    "HIMMELSKAMP":      {"type":  bool, "pos": 0}, # Include the Himmelskamp Stall delay? 
    "TOWERSHADOW":      {"type":  bool, "pos": 0}, # Include the tower shadow effect?
    "TOWERDRAG":        {"type": float, "pos": 0}, # Tower drag coefficient, [-]
    "WAKETYPE":         {"type":   int, "pos": 0}, # Wake type (0=free vortex wake, 1=unsteady BEM)
    "WAKEINTTYPE":      {"type":   int, "pos": 0}, # Wake integration type (0=EF, 1=ET, 2=PC, 3=PC2B)
    "WAKEROLLUP":       {"type":  bool, "pos": 0}, # Calculate wake self-induction?
    "TRAILINGVORT":     {"type":  bool, "pos": 0}, # Include trailing vortex elements?
    "SHEDVORT":         {"type":  bool, "pos": 0}, # Include shed vortex elements?
    "CONVECTIONTYPE":   {"type":   int, "pos": 0}, # Wake convection type (0=BL, 1=HH, 2=LOC)
    "WAKERELAXATION":   {"type": float, "pos": 0}, # Wake relaxation factor, [-]
    "FIRSTWAKEROW":     {"type": float, "pos": 0}, # First wake row length, [-]
    "MAXWAKESIZE":      {"type":   int, "pos": 0}, # Maximum number of wake elements, [-]
    "MAXWAKEDIST":      {"type": float, "pos": 0}, # Maximum wake distance from the rotor plane, [-]
    "WAKEREDUCTION":    {"type": float, "pos": 0}, # Wake reduction factor, [-]
    "WAKELENGTHTYPE":   {"type":   int, "pos": 0}, # Wake length type (0=counted in rotor revolutions, 1=counted in time steps)
    "CONVERSIONLENGTH": {"type": float, "pos": 0}, # Wake conversion length, [-]
    "NEARWAKELENGTH":   {"type": float, "pos": 0}, # Near wake length, [-]
    "ZONE1LENGTH":      {"type": float, "pos": 0}, # Wake zone 1 length, [-]
    "ZONE2LENGTH":      {"type": float, "pos": 0}, # Wake zone 2 length, [-]
    "ZONE3LENGTH":      {"type": float, "pos": 0}, # Wake zone 3 length, [-]
    "ZONE1FACTOR":      {"type":   int, "pos": 0}, # Wake zone 1 streamwise factor, [-]
    "ZONE2FACTOR":      {"type":   int, "pos": 0}, # Wake zone 2 streamwise factor, [-]
    "ZONE3FACTOR":      {"type":   int, "pos": 0}, # Wake zone 3 streamwise factor, [-]
    "ZONE1FACTOR_S":    {"type":   int, "pos": 0}, # Wake zone 1 spanwise factor, [-]
    "ZONE2FACTOR_S":    {"type":   int, "pos": 0}, # Wake zone 2 spanwise factor, [-]
    "ZONE3FACTOR_S":    {"type":   int, "pos": 0}, # Wake zone 3 spanwise factor, [-]
    "BOUNDCORERADIUS":  {"type": float, "pos": 0}, # Fixed core radius of the bound blade vortex, [-]
    "WAKECORERADIUS":   {"type": float, "pos": 0}, # Initial core radius of the free wake vortex, [-]
    "VORTEXVISCOSITY":  {"type": float, "pos": 0}, # Turbulent vortex viscosity, [-]
    "VORTEXSTRAIN":     {"type":  bool, "pos": 0}, # Calculate vortex strain?
    "MAXSTRAIN":        {"type":   int, "pos": 0}, # Maximum element strain, [-]
    "GAMMARELAXATION":  {"type": float, "pos": 0}, # Relaxation factor used in the gamma iteration, [-]
    "GAMMAEPSILON":     {"type": float, "pos": 0}, # Relative gamma convergence criteria, [-]
    "GAMMAITERATIONS":  {"type":   int, "pos": 0}, # Maximum number of gamma iterations, [-]
    "POLARDISC":        {"type":   int, "pos": 0}, # Polar discretization for the unsteady BEM, [-]
    "BEMTIPLOSS":       {"type":  bool, "pos": 0}, # Use BEM tip loss factor?
    "BEMSPEEDUP":       {"type": float, "pos": 0}, # Initial BEM convergence acceleration time, [s]
    "STRUCTURALFILE":   {"type":   str, "pos": 0}, # Input file for the structural model (leave blank if unused)
    "GEOMSTIFFNESS":    {"type":  bool, "pos": 0}, # Enable geometric stiffness?
    "AEROPANELLOADS":   {"type":  bool, "pos": 0}, # Enable distributed aero panel loads and gradients?
    "CONTROLLERTYPE":   {"type":   int, "pos": 0}, # Type of turbine controller (0=none, 1=BLADED, 2=DTU, 3=TUB)
    "CONTROLLERFILE":   {"type":   str, "pos": 0}, # Controller file name (leave blank if unused)
    "PARAMETERFILE":    {"type":   str, "pos": 0}, # Controller parameter file name (leave blank if unused)
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
        blade: Blade -- Blade object
        control : Control -- Control object
        path : str -- path to the .trb file
        structure : Structure -- Structure object
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

        # Check if file exists
        if os.path.isfile(path):
            self.path = path
        else:
            print(f"No file found at {path}!")
            sys.exit(1)

        # Read file
        self.read()

        # Add Aero, Structure, and Control objects
        self.blade = Blade(self.attributes["BLADEFILE"])
        self.structure = None
        self.control = None

    def read(self):
        """
        Read the .trb file.

        Arguments:
            None

        Returns:
            None
        """

        # Open file
        f = open(self.path, "r", encoding="utf-8")

        # Read attributes
        for key, value in TURBINE_DICT.items():
            self.attributes[key] = read(f, key, value["type"], value["pos"])

        self.attributes["BLADEFILE"]      = self.attributes["BLADEFILE"].replace("/", "\\")
        self.attributes["STRUCTURALFILE"] = self.attributes["STRUCTURALFILE"].replace("/", "\\")
        self.attributes["CONTROLLERFILE"] = self.attributes["CONTROLLERFILE"].replace("/", "\\")
        self.attributes["PARAMETERFILE"]  = self.attributes["PARAMETERFILE"].replace("/", "\\")

        self.attributes["BLADEFILE"]      = os.path.join(os.path.dirname(self.path), self.attributes["BLADEFILE"])
        self.attributes["STRUCTURALFILE"] = os.path.join(os.path.dirname(self.path), self.attributes["STRUCTURALFILE"])
        self.attributes["CONTROLLERFILE"] = os.path.join(os.path.dirname(self.path), self.attributes["CONTROLLERFILE"])
        self.attributes["PARAMETERFILE"]  = os.path.join(os.path.dirname(self.path), self.attributes["PARAMETERFILE"])

        self.attributes["SHAFTTILT"] = np.radians(self.attributes["SHAFTTILT"])
        self.attributes["ROTORCONE"] = np.radians(self.attributes["ROTORCONE"])
        self.attributes["XTILT"]     = np.radians(self.attributes["XTILT"])
        self.attributes["YTILT"]     = np.radians(self.attributes["YTILT"])

        # Close file
        f.close()

    def write(self):
        pass

if __name__ == "__main__":

    # Create a Turbine instance
    turbine = Turbine("data\\turbines\\DTU_10MW\\DTU_10MW_RWT.trb")

    # Print attributes
    for key, value in turbine.attributes.items():
        print(f"{key}: {value}")
