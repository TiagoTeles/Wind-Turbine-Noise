"""
Author:   T. Moreira da Fonte Fonseca Teles
Email:    tmoreiradafont@tudelft.nl
Date:     2025-11-10
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
        read_results -- read the simulation results
        get_results -- get the simulation results

    Attributes:
        path : str -- path to the .trb file
        n_blades : int -- number of blades, [-]
        n_panels : int -- number of panels, [-]
        rotor_overhang : float -- rotor overhang, [m]
        shaft_tilt : float -- shaft tilt angle, [rad]
        rotor_cone : float -- rotor cone angle, [rad]
        tower_height : float -- tower height, [m]
        blade: Blade -- blade object
        azimuth : np.ndarray -- azimuth angle, [rad]
        tip_speed_ratio : np.ndarray -- tip speed ratio, [-]
        yaw : np.ndarray -- yaw angle, [rad]
        power_coefficient : np.ndarray -- power coefficient, [-]
        torque_coefficient : np.ndarray -- torque coefficient, [-]
        thrust_coefficient : np.ndarray -- thrust coefficient, [-]
        pitch : np.ndarray -- pitch angle, [rad]
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

        # Read the file
        f = open(self.path, "r", encoding="utf-8")
        lines = f.readlines()

        # Set the turbine geometry
        self.n_blades = int(lines[12].split()[0])
        self.n_panels = int(lines[17].split()[0])
        self.rotor_overhang = float(lines[21].split()[0])
        self.shaft_tilt = float(lines[22].split()[0])
        self.rotor_cone = float(lines[23].split()[0])
        self.tower_height = float(lines[27].split()[0])

        # Convert the angles from [deg] to [rad]
        self.shaft_tilt = np.radians(self.shaft_tilt)
        self.rotor_cone = np.radians(self.rotor_cone)

        # Add the Blade object
        blade_directory = os.path.dirname(self.path)
        blade_name = lines[10].split()[0]
        blade_path = os.path.normpath(os.path.join(blade_directory, blade_name))
        self.blade = Blade(blade_path)

        # Close the file
        f.close()

        # Initialise the results attributes
        self.azimuth = None
        self.tip_speed_ratio = None
        self.yaw = None
        self.power_coefficient = None
        self.torque_coefficient  = None
        self.thrust_coefficient  = None
        self.pitch = None

    def read_results(self, results):
        """
        Read the simulation results.

        Parameters:
            results : pd.DataFrame -- simulation results

        Returns:
            None
        """

        # Read the Turbine results
        self.azimuth = results["Azimuthal~Angle~BLD_1~[deg]"].to_numpy()
        self.tip_speed_ratio = results["Tip~Speed~Ratio~[-]"].to_numpy()
        self.yaw = results["Yaw~Angle~[deg]"].to_numpy()
        self.power_coefficient = results["Power~Coefficient~[-]"].to_numpy()
        self.torque_coefficient = results["Torque~Coefficient~[-]"].to_numpy()
        self.thrust_coefficient = results["Thrust~Coefficient~[-]"].to_numpy()
        self.pitch = results["Pitch~Angle~BLD_1~[deg]"].to_numpy()

        # Convert the angles from [deg] to [rad]
        self.azimuth = np.radians(self.azimuth)
        self.yaw = np.radians(self.yaw)
        self.pitch = np.radians(self.pitch)

        # Read the Blade results
        self.blade.read_results(results)

    def get_results(self, key, azimuth):
        """
        Get the simulation results.

        Parameters:
            key : str -- member name
            azimuth : float -- azimuth angle, [rad]
        """

        # Check if the key is valid
        if key in ["azimuth", "tip_speed_ratio", "yaw", "power_coefficient", \
                   "torque_coefficient", "thrust_coefficient", "pitch"]:

            # Interpolate the results
            value = np.interp(azimuth, self.azimuth, getattr(self, key), period=2*np.pi)

        else:
            print("Invalid key!")
            sys.exit(1)

        return value
