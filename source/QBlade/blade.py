"""
Author:   T. Moreira da Fonte Fonseca Teles
Email:    tmoreiradafont@tudelft.nl
Date:     2025-09-22
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

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd

from source.QBlade.polar import Polar
from source.settings import SIMULATION_PATH, AR


class Blade():
    """
    A class to store the blade data.

    Methods:
        __init__ -- initialise the Blade class
        get -- get the blade property at a given radius
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

    def get(self, key, radius):
        """
        Get the blade geometry at a given radius.

        Parameters:
            key : str -- property key
            radius : np.ndarray -- radius, [m]

        Returns:
            value : np.ndarray -- property value
        """

        # Determine the value
        value = np.interp(radius, self.radius, getattr(self, key))

        return value

    def discretise(self, AR):
        """
        Discretise the blade into panels.

        Parameters:
            AR : float -- aspect ratio, [-]

        Returns:
            radius_p : np.ndarray -- panel radiuses, [m]
            span_p : np.ndarray -- panel spans, [m]
            chord_p : np.ndarray -- panel chords, [m]
        """

        radius = []
        chord = []

        # Start at the rotor blade tip
        radius.append(self.radius[-1])
        chord.append(self.chord[-1])

        # Loop until the blade root is reached
        while radius[-1] > self.radius[0]:

            # Assume an initial guess
            c_panel = chord[-1]
            AR_panel = 0.0

            # Loop until the AR is reached
            while np.abs(AR_panel - AR) > 1.0E-3:

                # Determine the radius and chord
                r_new = radius[-1] - AR * c_panel
                c_new = self.get("chord", r_new)

                # Determine the aspect ratio
                b_panel = radius[-1] - r_new
                c_panel = (chord[-1] + c_new) / 2.0
                AR_panel = b_panel / c_panel

            # Save the new radius and chord
            radius.append(r_new)
            chord.append(c_new)

        # Include the rotor blade root
        radius[-1] = self.radius[0]
        chord[-1] = self.chord[0]

        # Invert the radius and chord lists
        radius = np.array(radius[::-1])
        chord = np.array(chord[::-1])

        # Determine the panel radius and chord
        radius_p = (radius[1:] + radius[:-1]) / 2.0
        chord_p = (chord[1:] + chord[:-1]) / 2.0

        # Determine the panel span
        span_p = np.diff(radius)

        return radius_p, span_p, chord_p

if __name__ == "__main__":

    # Import the Simulation class
    from source.QBlade.simulation import Simulation

    # Create the Blade object
    simulation = Simulation(SIMULATION_PATH)
    turbine = simulation.turbine
    blade = turbine.blade

    # Discretise the blade
    radius, span, chord = blade.discretise(AR)

    # Show the blade discretisation
    plt.bar(radius, chord, width=span, color="tab:blue", edgecolor="black", label="Panels", zorder=2)
    plt.plot(blade.radius, blade.chord, color="tab:orange", label="Chord", zorder=3)

    plt.xlabel("Radius, [m]")
    plt.ylabel("Chord, [m]")

    plt.xlim(0.0, 146.2)
    plt.ylim(0.0, 8.0)

    plt.legend()
    plt.grid(zorder=0)
    plt.show()

    # Determine the number of panels for different ARs
    AR_list = []
    n_list = []

    for AR in np.linspace(1.0, 10.0, 1000):

        # Discretise the blade
        radius, span, chord = blade.discretise(AR)

        # Save the AR and number of panels
        AR_list.append(AR)
        n_list.append(len(radius))

    # Show the number of panels
    plt.plot(AR_list, n_list)

    plt.xlabel("Panel Aspect Ratio, [-]")
    plt.ylabel("Number of Panels, [-]")

    plt.xlim(1.0, 10.0)
    plt.ylim(0.0, 40.0)

    plt.grid()
    plt.show()
