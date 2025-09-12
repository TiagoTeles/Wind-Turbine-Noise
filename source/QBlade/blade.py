"""
Author:   T. Moreira da Fonte Fonseca Teles
Email:    tmoreiradafont@tudelft.nl
Date:     2025-09-10
License:  GNU GPL 3.0

Store the blade data.

Classes:
    Blade

Functions:
    None

Exceptions:
    None
"""

import glob
import os
import sys

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
        read -- read the .bld file
        discretise -- discretise the blade into panels

    Attributes:
        path : str -- path to the .bld file
        geometry : pd.DataFrame -- blade geometry
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

        # Check if the file exists
        if not os.path.isfile(path):
            print(f"No file found at {path}!")
            sys.exit(1)

        # Read the file
        self.read()

    def read(self):
        """
        Read the .bld file.

        Parameters:
            None

        Returns:
            None
        """

        # Read the blade geometry
        self.geometry = pd.read_csv(self.path, delimiter=r"\s+", names=["radius", "chord", "twist", \
                                    "offset_x", "offset_y", "pitch_axis", "polar_path"], skiprows=16)

        # Convert the twist from [deg] to [rad]
        self.geometry["twist"] = np.radians(self.geometry["twist"])

        # Add the Polar object
        self.geometry["polar"] = None

        for index, row in self.geometry.iterrows():
            polar_path = os.path.join(os.path.dirname(self.path), row["polar_path"])
            self.geometry.at[index, "polar"] = Polar(polar_path)

    def discretise(self, AR):
        """
        Discretise the blade into panels.

        Parameters:
            AR : float -- aspect ratio of the panels, [-]

        Returns:
            radius_p : list -- list of panel radiuses, [m]
            span_p : list -- list of panel spans, [m]
            chord_p : list -- list of panel chords, [m]
        """

        radius = []
        chord = []

        # Start the discretisation at the tip
        radius.append(self.geometry["radius"].iat[-1])
        chord.append(self.geometry["chord"].iat[-1])

        # Iterate though the entire blade
        while True:

            # Set the outer panel radius and chord
            radius_0 = radius[-1]
            chord_0 = chord[-1]

            # Guess the inner panel radius and chord
            radius_1 = radius_0
            chord_1 = chord_0

            # Guess the MAC
            MAC_01 = (chord_0 + chord_1) / 2.0

            # Iterate until the AR is correct
            while True:

                # Update the inner panel radius and chord
                radius_1 = radius_0 - AR * MAC_01
                chord_1 = np.interp(radius_1, self.geometry["radius"], self.geometry["chord"])

                # Update the MAC
                MAC_01 = (chord_0 + chord_1) / 2.0

                # Determine the current AR
                AR_01 = (radius_0 - radius_1) / MAC_01

                # Check if the AR is correct
                if np.abs(AR_01 - AR) < 1.0E-6:

                    # Set the inner panel radius and chord
                    radius.append(radius_1)
                    chord.append(chord_1)

                    break

            # Check if the end of the blade is reached
            if radius[-1] < self.geometry.at[0, "radius"]:

                # Set the last panel radius and chord
                radius[-1] = self.geometry.at[0, "radius"]
                chord[-1] = self.geometry.at[0, "chord"]

                break

        # Invert the radius and chord lists
        radius = np.array(radius[::-1])
        chord =  np.array(chord[::-1])

        # Determine the panel radius and chord
        radius_p = (radius[1:] + radius[:-1]) / 2.0
        chord_p = (chord[1:] + chord[:-1]) / 2.0

        # Determine the panel span
        span_p = np.diff(radius)

        return radius_p, span_p, chord_p

if __name__ == "__main__":

    # Determine the blade directory
    blade_path = glob.glob(os.path.dirname(SIMULATION_PATH) + "\\**\\*.bld", recursive=True)[0]

    # Create the Blade object
    blade = Blade(blade_path)

    # Discretise the blade
    radius, span, chord = blade.discretise(AR)

    # Show the blade discretisation
    plt.bar(radius, chord, width=span, color="tab:blue", edgecolor="black", label="Panels", zorder=2)
    plt.plot(blade.geometry["radius"], blade.geometry["chord"], color="tab:orange", label="Chord", zorder=3)

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

    plt.xlabel("Aspect Ratio, [-]")
    plt.ylabel("Number of Panels, [-]")

    plt.xlim(1.0, 10.0)
    plt.ylim(0.0, 40.0)

    plt.grid()
    plt.show()
