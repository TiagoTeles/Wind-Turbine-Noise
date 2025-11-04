"""
Author:   T. Moreira da Fonte Fonseca Teles
Email:    tmoreiradafont@tudelft.nl
Date:     2025-11-04
License:  GNU GPL 3.0

Store the blade data.

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
import scipy as sp

from source.QBlade.polar import Polar


class Blade():
    """
    A class to store the blade data.

    Methods:
        __init__ -- initialise the Blade class
        get_geometry -- get the blade geometry
        read_results -- read the simulation results
        get_results -- get the simulation results
        discretise -- discretise the blade into panels

    Attributes:
        path : str -- path to the .bld file
        radius : np.ndarray -- spanwise position, [m]
        chord : np.ndarray -- chord length, [m]
        twist : np.ndarray -- twist angle, [rad]
        offset_x : np.ndarray -- in-plane offset, [m]
        offset_y : np.ndarray -- out-of-plane offset, [m]
        pitch_axis : np.ndarray -- pitch axis position, [-]
        polar : np.ndarray -- polar data
        thickness_01 : np.ndarray -- airfoil thickness at x/c = 0.01 [-]
        thickness_10 : np.ndarray -- airfoil thickness at x/c = 0.10 [-]
        azimuth : np.ndarray -- azimuth angle, [rad]
        angle_of_attack : np.ndarray -- angle of attack, [rad]
        panel_radius : np.ndarray -- panel radius, [m]
        axial_force : np.ndarray -- axial force per unit span, [N/m]
        tangential_force : np.ndarray -- tangential force per unit span, [N/m]
        velocity : np.ndarray -- velocity, [m/s]
        beam_radius : np.ndarray -- beam radius, [m]
        axial_deflection : np.ndarray -- axial deflection, [m]
        radial_twist : np.ndarray -- radial twist, [rad]
    """

    def __init__(self, path):
        """
        Initialise the Blade class.

        Parameters:
            path : str -- path to the .bld file

        Returns:
            None
        """

        self.path = path

        # Read the blade geometry
        geometry = pd.read_csv(self.path, delimiter=r"\s+", names=["radius", "chord", "twist", \
                               "offset_x", "offset_y", "pitch_axis", "polar_path"], skiprows=16)

        # Add the blade geometry
        self.radius = geometry["radius"].to_numpy()
        self.chord = geometry["chord"].to_numpy()
        self.twist = geometry["twist"].to_numpy()
        self.offset_x = geometry["offset_x"].to_numpy()
        self.offset_y = geometry["offset_y"].to_numpy()
        self.pitch_axis = geometry["pitch_axis"].to_numpy()

        # Convert the twist from [deg] to [rad]
        self.twist = np.radians(self.twist)

        # Add the Polar object
        self.polar = np.empty(self.radius.shape, dtype=Polar)

        for i in range(len(self.radius)):
            polar_path = os.path.normpath(os.path.join(os.path.dirname(self.path), geometry["polar_path"].iat[i]))
            self.polar[i] = Polar(polar_path)

        # Add the airfoil thicknesses
        self.thickness_01 = np.empty(self.radius.shape)
        self.thickness_10 = np.empty(self.radius.shape)

        for i in range(len(self.radius)):
            self.thickness_01[i] = self.polar[i].airfoil.thickness(0.01)
            self.thickness_10[i] = self.polar[i].airfoil.thickness(0.10)

        # Initialise the results attributes
        self.azimuth = None
        self.angle_of_attack = None
        self.panel_radius = None
        self.axial_force = None
        self.tangential_force = None
        self.velocity = None
        self.beam_radius = None
        self.axial_deflection = None
        self.radial_twist = None

    def get_geometry(self, key, radius):
        """
        Get the blade geometry.

        Parameters:
            key : str -- member name
            radius : np.ndarray -- radius, [m]

        Returns:
            value : np.ndarray -- property value
        """

        # Determine the value
        value = np.interp(radius, self.radius, getattr(self, key))

        return value

    def read_results(self, results):
        """
        Read the simulation results.

        Parameters:
            results : pd.DataFrame -- simulation results

        Returns:
            None
        """

        # Determine the number of timesteps
        n_timesteps = len(results["Time~[s]"])

        # Determine the number of panels and beams
        n_panels = 0
        n_beams = 0

        for column in results.columns.values:
            if "Radius~BLD_1~PAN_" in column:
                n_panels += 1

            if "X_l~For.~BLD_1~pos~" in column:
                n_beams += 1

        # Read the panel radius and beam radius
        self.panel_radius = np.empty(n_panels)
        self.beam_radius = np.linspace(self.radius[0], self.radius[-1], n_beams, endpoint=True)

        for i in range(n_panels):
            self.panel_radius[i] = results[f"Radius~BLD_1~PAN_{i}~[m]"].iat[0]

        # Read the azimuth angle
        self.azimuth = results["Azimuthal~Angle~BLD_1~[deg]"].to_numpy()

        # Read the aerodynamic blade distributions
        self.angle_of_attack = np.empty((n_timesteps, n_panels))
        self.axial_force = np.empty((n_timesteps, n_panels))
        self.tangential_force = np.empty((n_timesteps, n_panels))
        self.velocity = np.empty((n_timesteps, n_panels))

        for i in range(n_panels):
            self.angle_of_attack[:, i] = results[f"Angle~of~Attack~at~0.25c~BLD_1~PAN_{i}~[deg]"]
            self.axial_force[:, i] = results[f"Normal~Force~BLD_1~PAN_{i}~[N/m]"]
            self.tangential_force[:, i] = results[f"Tangential~Force~BLD_1~PAN_{i}~[N/m]"]
            self.velocity[:, i] = results[f"Total~Velocity~BLD_1~PAN_{i}~[m/s]"]

        # Read the structural blade distributions
        self.axial_deflection = np.empty((n_timesteps, n_beams))
        self.radial_twist = np.empty((n_timesteps, n_beams))

        for i in range(n_beams):
            self.axial_deflection[:, i] = results[f"X_b~Trl.Def.~BLD_1~pos~{(i / (n_beams - 1)):.3f}~[m]"]
            self.radial_twist[:, i] = results[f"Z_b~Rot.Def.~BLD_1~pos~{(i / (n_beams - 1)):.3f}~[deg]"]

        # Convert the angles from [deg] to [rad]
        self.azimuth = np.radians(self.azimuth)
        self.angle_of_attack = np.radians(self.angle_of_attack)
        self.radial_twist = np.radians(self.radial_twist)

    def get_results(self, key, azimuth, radius):
        """
        Get the simulation results.

        Parameters:
            key : str -- member key
            azimuth : np.array -- azimuth angle, [rad]
            radius : np.ndarray -- radius, [m]
        """
    
        # Clamp the azimuth angle
        azimuth = azimuth % (2 * np.pi)

        # Check if the key is valid
        if key in ["angle_of_attack", "axial_force", "tangential_force", "velocity"]:

            # Create a 2D interpolator
            interpolator = sp.interpolate.RegularGridInterpolator((self.azimuth, self.panel_radius), \
                           getattr(self, key), bounds_error=False, fill_value=None)

        elif key in ["axial_deflection", "radial_twist"]:

            # Create a 2D interpolator
            interpolator = sp.interpolate.RegularGridInterpolator((self.azimuth, self.beam_radius), \
                           getattr(self, key), bounds_error=False, fill_value=None)

        else:
            print("Invalid key!")
            sys.exit(1)

        # Interpolate the results
        value = interpolator((azimuth, radius))

        return value
