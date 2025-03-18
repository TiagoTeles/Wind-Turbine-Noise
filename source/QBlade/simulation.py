"""
Author:   T. Moreira da Fonte Fonseca Teles
Email:    tmoreiradafont@tudelft.nl
Date:     2025-03-18
License:  GNU GPL 3.0

Store the data from .sim files.

Classes:
    Simulation

Functions:
    None

Exceptions:
    None
"""

import os
import sys

import numpy as np

from QBlade.misc import read, write
from QBlade.turbine import Turbine


SIMULATION_DICT = {
    "OBJECTNAME":      {"type":   str}, # Name of the simulation object
    "ISOFFSHORE":      {"type":   int}, # Use a number (0=onshore, 1=offshore)
    "TURBFILE":        {"type":   str}, # Turbine definition file path
    "TURBNAME":        {"type":   str}, # Name of the turbine in the simulation
    "INITIAL_YAW":     {"type": float}, # Initial turbine yaw, [rad]
    "INITIAL_PITCH":   {"type": float}, # Initial collective blade pitch, [rad]
    "INITIAL_AZIMUTH": {"type": float}, # Initial azimuthal rotor angle, [rad]
    "STRSUBSTEP":      {"type":   int}, # Number of structural substeps per timestep, [-]
    "RELAXSTEPS":      {"type":   int}, # Number of initial static structural relaxation steps, [-]
    "RELAXTIMESTEP":   {"type": float}, # Timestep for the initial static structural relaxation steps, [s]
    "PRESCRIBETYPE":   {"type":   int}, # Rotor RPM prescribe type (0=ramp-up; 1=whole sim; 2=no RPM prescibed) 
    "RPMPRESCRIBED":   {"type": float}, # Prescribed rotor RPM, [min^-1]
    "STRITERATIONS":   {"type":   int}, # Number of iterations for the time integration, [-]
    "MODNEWTONITER":   {"type":   int}, # Use the modified newton iteration?
    "INCLUDEAERO":     {"type":   int}, # Include aerodynamic forces?
    "INCLUDEHYDRO":    {"type":   int}, # Include hydrodynamic forces?
    "GLOBPOS_X":       {"type": float}, # Global x-position of the turbine, [m]
    "GLOBPOS_Y":       {"type": float}, # Global y-position of the turbine, [m]
    "GLOBPOS_Z":       {"type": float}, # Global z-position of the turbine, [m]
    "GLOBROT_X":       {"type": float}, # Global x-rotation of the turbine, [rad]
    "GLOBROT_Y":       {"type": float}, # Global y-rotation of the turbine, [rad]
    "GLOBROT_Z":       {"type": float}, # Global z-rotation of the turbine, [rad]
    "EVENTFILE":       {"type":   str}, # File containing fault event definitions (leave blank if unused)
    "LOADINGFILE":     {"type":   str}, # Loading file name (leave blank if unused)
    "SIMFILE":         {"type":   str}, # Simulation file name (leave blank if unused)
    "MOTIONFILE":      {"type":   str}, # Prescribed motion file name (leave blank if unused)
    "TIMESTEP":        {"type": float}, # Timestep size, [s]
    "NUMTIMESTEPS":    {"type":   int}, # Number of timesteps, [-]
    "RAMPUP":          {"type": float}, # Rampup time for the structural model, [s]
    "ADDDAMP":         {"type": float}, # Initial time with additional damping, [s]
    "ADDDAMPFACTOR":   {"type": float}, # Factor used to increase the damping of all components, [-]
    "WAKEINTERACTION": {"type": float}, # Wake interaction start time for multi-turbine simulations, [s]
    "WNDTYPE":         {"type":   int}, # Use a number (0=steady, 1=windfield, 2=hubheight)
    "WNDNAME":         {"type":   str}, # Filename of the turbsim input file (leave blank if unused)
    "STITCHINGTYPE":   {"type":   int}, # Windfield stitching type (0=periodic; 1=mirror)
    "WINDAUTOSHIFT":   {"type":  bool}, # Windfield shifting automatically based on rotor diameter?
    "SHIFTTIME":       {"type": float}, # Windfield shift time, [s]
    "MEANINF":         {"type": float}, # Mean inflow velocity, [m/s]
    "HORANGLE":        {"type": float}, # Horizontal inflow angle, [rad]
    "VERTANGLE":       {"type": float}, # Vertical inflow angle, [rad]
    "PROFILETYPE":     {"type":   int}, # Type of wind profile used (0=Power Law; 1=Logarithmic)
    "SHEAREXP":        {"type": float}, # Shear exponent if using a power law profile, [-]
    "ROUGHLENGTH":     {"type": float}, # Roughness length if using a log profile, [m]
    "DIRSHEAR":        {"type": float}, # A value for the directional shear, [rad/m]
    "REFHEIGHT":       {"type": float}, # Reference height used to construct the BL profile, [m]
    "WATERDEPTH":      {"type": float}, # Water depth, [m]
    "WAVEFILE":        {"type":   str}, # Path to the wave file (leave blank if unused)
    "WAVESTRETCHING":  {"type":   int}, # Type of wave stretching (0=vertical, 1=wheeler, 2=extrapolation, 3=none)
    "SEABEDSTIFF":     {"type": float}, # Vertical seabed stiffness, [N/m^3]
    "SEABEDDAMP":      {"type": float}, # Damping factor for the vertical seabed stiffness evaluation, [-]
    "SEABEDSHEAR":     {"type": float}, # Factor for the evaluation of shear forces, [-]
    "SURF_CURR_U":     {"type": float}, # Near surface current velocity, [m/s]
    "SURF_CURR_DIR":   {"type": float}, # Near surface current direction, [rad]
    "SURF_CURR_DEPTH": {"type": float}, # Near surface current depth, [m]
    "SUB_CURR_U":      {"type": float}, # Sub surface current velocity, [m/s]
    "SUB_CURR_DIR":    {"type": float}, # Sub surface current direction, [rad]
    "SUB_CURR_EXP":    {"type": float}, # Sub surface current exponent, [-]
    "SHORE_CURR_U":    {"type": float}, # Near shore current velocity, [m/s]
    "SHORE_CURR_DIR":  {"type": float}, # Near shore current direction, [rad]
    "MOORINGSYSTEM":   {"type":   str}, # Path to the global mooring system file (leave blank if unused)
    "DWMSUMTYPE":      {"type":   int}, # Dynamic wake meandering wake summation type (0=DOMINANT; 1=QUADRATIC; 2=LINEAR)
    "DENSITYAIR":      {"type": float}, # Air density, [kg/m^3]
    "VISCOSITYAIR":    {"type": float}, # Air kinematic viscosity, [m^2/s]
    "DENSITYWATER":    {"type": float}, # Water density, [kg/m^3]
    "VISCOSITYWATER":  {"type": float}, # Water kinematic viscosity, [m^2/s]
    "GRAVITY":         {"type": float}, # Gravity constant, [m/s^2]
    "STOREFROM":       {"type": float}, # Simulation stores data from this point in time, [s]
    "STOREREPLAY":     {"type":  bool}, # Store a replay of the simulation?
    "STOREAERO":       {"type":  bool}, # Store the aerodynamic data?
    "STOREBLADE":      {"type":  bool}, # Store the local aerodynamic blade data?
    "STORESTRUCT":     {"type":  bool}, # Store the structural data?
    "STORESIM":        {"type":  bool}, # Store the simulation (performance) data?
    "STOREHYDRO":      {"type":  bool}, # Store the controller data?
    "STORECONTROLLER": {"type":  bool}, # Store the controller data?
    "STOREDWM":        {"type":  bool}, # Store the dynamic wake meandering data?
    "CALCMODAL":       {"type":  bool}, # Perform a modal analysis?
    "USEMBC":          {"type":  bool}, # Apply the multi blade coordinate transformation during the modal analysis?
    "MINFREQ":         {"type": float}, # Store Eigenvalues, starting with this frequency, [Hz]
    "DELTAFREQ":       {"type": float}, # Omit Eigenvalues that are closer spaced than this value, [Hz]
    "NUMFREQ":         {"type": float}, # Set the number of Eigenmodes and Eigenvalues that will be stored, [-]
    }

class Simulation:
    """
    A class to store the simulation data.

    Methods:
        __init__ -- initialise the Simulation class
        read -- read the .sim file
        write -- write to the .sim file

    Attributes:
        attributes : dict -- dictionary of attributes
        path : str -- path to the .trb file
        turbine : Turbine -- turbine object
    """

    def __init__(self, path):
        """
        Initialise the Simulation class.

        Arguments:
            path : str -- path to the .sim file

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

        # Add the Turbine object
        turbine_path = os.path.join(os.path.dirname(self.path), self.attributes["TURBFILE"])
        self.turbine = Turbine(turbine_path)

    def read(self):
        """
        Read the .sim file.

        Arguments:
            None

        Returns:
            None
        """

        # Open the file
        f = open(self.path, "r", encoding="utf-8")

        # Read the attributes
        for key, value in SIMULATION_DICT.items():
            self.attributes[key] = read(f, key, value["type"])

        # Format the attributes
        self.attributes["TURBFILE"]      = self.attributes["TURBFILE"].replace("/", "\\")
        self.attributes["EVENTFILE"]     = self.attributes["EVENTFILE"].replace("/", "\\")
        self.attributes["LOADINGFILE"]   = self.attributes["LOADINGFILE"].replace("/", "\\")
        self.attributes["SIMFILE"]       = self.attributes["SIMFILE"].replace("/", "\\")
        self.attributes["MOTIONFILE"]    = self.attributes["MOTIONFILE"].replace("/", "\\")
        self.attributes["WNDNAME"]       = self.attributes["WNDNAME"].replace("/", "\\")
        self.attributes["WAVEFILE"]      = self.attributes["WAVEFILE"].replace("/", "\\")
        self.attributes["MOORINGSYSTEM"] = self.attributes["MOORINGSYSTEM"].replace("/", "\\")

        self.attributes["INITIAL_YAW"]     = np.radians(self.attributes["INITIAL_YAW"])
        self.attributes["INITIAL_PITCH"]   = np.radians(self.attributes["INITIAL_PITCH"])
        self.attributes["INITIAL_AZIMUTH"] = np.radians(self.attributes["INITIAL_AZIMUTH"])
        self.attributes["GLOBROT_X"]       = np.radians(self.attributes["GLOBROT_X"])
        self.attributes["GLOBROT_Y"]       = np.radians(self.attributes["GLOBROT_Y"])
        self.attributes["GLOBROT_Z"]       = np.radians(self.attributes["GLOBROT_Z"])
        self.attributes["DIRSHEAR"]        = np.radians(self.attributes["DIRSHEAR"])
        self.attributes["SURF_CURR_DIR"]   = np.radians(self.attributes["SURF_CURR_DIR"])
        self.attributes["SUB_CURR_DIR"]    = np.radians(self.attributes["SUB_CURR_DIR"])
        self.attributes["SHORE_CURR_DIR"]  = np.radians(self.attributes["SHORE_CURR_DIR"])
        self.attributes["HORANGLE"]        = np.radians(self.attributes["HORANGLE"])
        self.attributes["VERTANGLE"]       = np.radians(self.attributes["VERTANGLE"])

        # Close the file
        f.close()

    def write(self, key, value):

        # Set the value in the attributes dictionary
        self.attributes[key] = value

        # Set the value in the .sim file
        f = open(self.path, "r+", encoding="utf-8")
        write(f, key, value)
        f.close()
