"""
Author:   T. Moreira da Fonte Fonseca Teles
Email:    tmoreiradafont@tudelft.nl
Date:     2025-03-21
License:  GNU GPL 3.0

Store the data from .sim files and run QBlade.

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
import pandas as pd

from QBlade.misc import read
from QBlade.qblade import QBlade
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
    A class to store the simulation data and run QBlade.

    Methods:
        __init__ -- initialise the Simulation class
        read -- read the .sim file
        initialise -- initialise the QBlade simulation
        run -- run the QBlade simulation
        close -- close the QBlade simulation

    Attributes:
        qblade : QBlade -- QBlade library
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
        self.cl_device = None
        self.dll_path = None
        self.group_size = None
        self.n_timestep = None
        self.qblade = None
        self.sim_path = path

        # Check if the file exists
        if not os.path.isfile(path):
            print(f"No file found at {path}!")
            sys.exit(1)

        # Read the file
        self.read()

        # Add the Turbine object
        turbine_path = os.path.join(os.path.dirname(self.sim_path), self.attributes["TURBFILE"])
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
        f = open(self.sim_path, "r", encoding="utf-8")

        # Read the attributes
        for key, value in SIMULATION_DICT.items():
            self.attributes[key] = read(f, key, value["type"])

        # Format the attributes
        self.attributes["TURBFILE"] = self.attributes["TURBFILE"].replace("/", "\\")
        self.attributes["EVENTFILE"] = self.attributes["EVENTFILE"].replace("/", "\\")
        self.attributes["LOADINGFILE"] = self.attributes["LOADINGFILE"].replace("/", "\\")
        self.attributes["SIMFILE"] = self.attributes["SIMFILE"].replace("/", "\\")
        self.attributes["MOTIONFILE"] = self.attributes["MOTIONFILE"].replace("/", "\\")
        self.attributes["WNDNAME"] = self.attributes["WNDNAME"].replace("/", "\\")
        self.attributes["WAVEFILE"] = self.attributes["WAVEFILE"].replace("/", "\\")
        self.attributes["MOORINGSYSTEM"] = self.attributes["MOORINGSYSTEM"].replace("/", "\\")

        self.attributes["INITIAL_YAW"] = np.radians(self.attributes["INITIAL_YAW"])
        self.attributes["INITIAL_PITCH"] = np.radians(self.attributes["INITIAL_PITCH"])
        self.attributes["INITIAL_AZIMUTH"] = np.radians(self.attributes["INITIAL_AZIMUTH"])
        self.attributes["GLOBROT_X"] = np.radians(self.attributes["GLOBROT_X"])
        self.attributes["GLOBROT_Y"] = np.radians(self.attributes["GLOBROT_Y"])
        self.attributes["GLOBROT_Z"] = np.radians(self.attributes["GLOBROT_Z"])
        self.attributes["DIRSHEAR"] = np.radians(self.attributes["DIRSHEAR"])
        self.attributes["SURF_CURR_DIR"] = np.radians(self.attributes["SURF_CURR_DIR"])
        self.attributes["SUB_CURR_DIR"] = np.radians(self.attributes["SUB_CURR_DIR"])
        self.attributes["SHORE_CURR_DIR"] = np.radians(self.attributes["SHORE_CURR_DIR"])
        self.attributes["HORANGLE"] = np.radians(self.attributes["HORANGLE"])
        self.attributes["VERTANGLE"] = np.radians(self.attributes["VERTANGLE"])

        # Close the file
        f.close()

    def initialise(self, path, cl_device, group_size):
        """
        Initialise the QBlade simulation.

        Arguments:
            path : str -- path to the QBlade.dll file
            cl_device : int -- OpenCL device
            group_size : int -- OpenCL group size

        Returns:
            None
        """

        self.dll_path = path
        self.cl_device = cl_device
        self.group_size = group_size

        # Load the QBlade library
        self.qblade = QBlade(path)
        self.qblade.createInstance(cl_device, group_size)

        # Setup the simulation
        self.qblade.loadSimDefinition(self.sim_path.encode("utf-8"))
        self.qblade.initializeSimulation()

    def run(self, n_timestep):
        """
        Run the QBlade simulation.

        Arguments:
            n_timestep : int -- number of timesteps

        Returns:
            None
        """

        self.n_timestep = n_timestep

        # Run the simulation
        for i in range(n_timestep):

            # Advance the simulation one timestep
            success = self.qblade.advanceTurbineSimulation()

            # Ensure the simulation step was successful
            if not success:
                print(f"Simulation failed at timestep {i}!")
                sys.exit(1)

        # Determine the results path
        results_dir = os.path.dirname(self.sim_path)
        results_name = os.path.splitext(os.path.basename(self.sim_path))[0]

        # Save the simulation results
        self.qblade.exportResults(0, results_dir.encode("utf-8"), results_name.encode("utf-8"), b"")

    def results(self, timestep, blade, c_0, cutoff, probe_top, probe_bot, it):
        """
        Read the simulation results.

        Arguments:
            timestep : int -- timestep index
            blade : int -- blade index
            c_0 : float -- speed of sound, [m]
            cutoff : float -- radial cutoff, [-]
            probe_top : np.array -- probe location at the top surface, [-]
            probe_bot : np.array -- probe location at the bottom surface, [-]
            it : int -- number of iterations for XFoil, [-]

        Returns:
            results : dict -- dictionary of simulation results
        """
        # Determine the results path
        results_dir = os.path.dirname(self.sim_path)
        results_name = os.path.splitext(os.path.basename(self.sim_path))[0]
        results_path = os.path.join(results_dir, results_name + ".txt")

        # Read the simulation results
        data = pd.read_csv(results_path, skiprows=2, delimiter="\t").iloc[timestep]
        results = {}

        # Read the turbine attributes
        pitch = np.radians(data[f"Pitch_Angle_Blade_{blade}_[deg]"])
        yaw = np.radians(data["Yaw_Angle_[deg]"])

        results["pitch"] = pitch
        results["yaw"] = yaw

        # Read the blade distributions
        n_panels = self.turbine.attributes["NUMPANELS"]

        aoa = np.zeros(n_panels)
        pos = np.zeros(n_panels)
        U = np.zeros(n_panels)
        Re = np.zeros(n_panels)

        for i in range(n_panels):
            aoa[i] = np.radians(data[f"Angle_of_Attack_at_0.25c_Blade_{blade}_PAN_{i}_[deg]"])
            pos[i] = data[f"Radius_Blade_{blade}_PAN_{i}_[m]"]
            U[i] = data[f"Total_Velocity_Blade_{blade}_PAN_{i}_[m/s]"]
            Re[i] = data[f"Reynolds_Number_Blade_{blade}_PAN_{i}_[-]"]

        results["aoa"] = aoa
        results["pos"] = pos
        results["U"] = U
        results["Re"] = Re

        # Interpolate the blade distributions
        radiuses = np.array(self.turbine.blade.data["pos"])
        chords = np.array(self.turbine.blade.data["chord"])
        twists    = np.array(self.turbine.blade.data["twist"])
        offsets_x = np.array(self.turbine.blade.data["offset_x"])
        offsets_y = np.array(self.turbine.blade.data["offset_y"])
        p_axes   = np.array(self.turbine.blade.data["p_axis"])

        chord = np.interp(pos, radiuses, chords)
        twist = np.interp(pos, radiuses, twists)
        offset_x = np.interp(pos, radiuses, offsets_x)
        offset_y = np.interp(pos, radiuses, offsets_y)
        p_axis = np.interp(pos, radiuses, p_axes)

        results["chord"] = chord
        results["twist"] = twist
        results["offset_x"] = offset_x
        results["offset_y"] = offset_y
        results["p_axis"] = p_axis

        # Determine the panel spanwise distributions
        if self.turbine.attributes["DISCTYPE"] == 0:
            spans = np.diff(radiuses)
            span = np.interp(pos, radiuses, spans)

        elif self.turbine.attributes["DISCTYPE"] == 1:
            span = np.ones(n_panels) * (radiuses[-1] - radiuses[0]) / n_panels

        elif self.turbine.attributes["DISCTYPE"] == 2:
            r_virtual = np.concat((np.array([radiuses[0]]), pos, np.array([radiuses[-1]])))
            span = (r_virtual[2:n_panels+2] - r_virtual[0:n_panels]) / 2

        else:
            print("Invalid DISCTYPE!")
            sys.exit(1)

        results["span"] = span

        # Determine the airfoil thickness distributions
        tc_01 = np.zeros(n_panels)
        tc_10 = np.zeros(n_panels)

        for i in range(n_panels):
            airfoil = self.turbine.blade.interpolate(pos[i])
            tc_01[i] = airfoil.thickness(0.01)
            tc_10[i] = airfoil.thickness(0.10)

        results["tc_01"] = tc_01
        results["tc_10"] = tc_10

        # Determine the airfoil displacement thickness distributions
        delta_star_top = np.zeros(n_panels)
        delta_star_bot = np.zeros(n_panels)

        for i in range(n_panels):
            airfoil = self.turbine.blade.interpolate(pos[i])

            delta_star_top[i], delta_star_bot[i] = \
                self.turbine.blade.displacement_thickness(pos[i], Re[i], U[i]/c_0, aoa[i], cutoff, probe_top, probe_bot, it)

        results["delta_star_top"] = delta_star_top
        results["delta_star_bot"] = delta_star_bot

        return results

    def close(self):
        """
        Close the QBlade simulation.

        Arguments:
            None

        Returns:
            None
        """

        # Unload the QBlade library
        self.qblade.unload_library()
        self.qblade = None
