""" 
Author:   T. Moreira da Fonte Fonseca Teles
Email:    tmoreiradafont@tudelft.nl
Date:     2025-02-14
License:  GNU GPL 3.0

Store aerodynamic data.

Classes:
    Aero

Functions:
    None

Exceptions:
    None
"""

import os
import sys
import numpy as np
import pandas as pd

from airfoil import Airfoil
from misc import parse
from polar import Polar


class Aero():
    """ 
    A class to store the aerodynamic data. 
    
    Methods:
        __init__ -- parse the blade file
        
    Attributes:
        airfoils : list -- a list of Airfoil objects
        airfoil_dir : str -- the path to the airfoil folder
        data : pd.DataFrame -- the blade data
        inverted_foils : bool -- invert the airfoils? (only VAWT)
        path : str -- the path to the .bld file
        polars : list -- a list of Polar objects
        polar_dir : str -- the path to the polar folder
        name : str -- the name of the blade object
        n_blades : int -- the number of blades
        type : str -- the rotor type
    """

    def __init__(self, path):
        """
        Initializes the Aero object.

        Methods:
            __init__ -- parse the blade file
            
        Attributes:
            path : str -- the path to the .bld file
        """

        self.path = path
        self.airfoils = []
        self.polars = []

        # Open file
        if os.path.isfile(path):
            f = open(path, "r", encoding="utf-8")
        else:
            print(f"No blade file found at {path}!")
            sys.exit(1)

        # Parse data in file
        self.name           = parse(f,    "OBJECTNAME", 0,  str)
        self.type           = parse(f,     "ROTORTYPE", 0,  str)
        self.inverted_foils = parse(f, "INVERTEDFOILS", 0, bool)
        self.n_blades       = parse(f,     "NUMBLADES", 0,  int)

        # Read pos, chord, twist, offset_x, offset_y, p_axis, and polar file
        self.data = pd.read_csv(f, names=['pos', 'chord', 'twist', 'offset_x', 'offset_y', \
                                          'p_axis', 'polar_path'], delimiter=r"\s+", skiprows=3)

        # Format parsed data
        self.data['twist'] = np.radians(self.data['twist'])
        self.data['polar_path'] = os.path.dirname(path) + os.sep + self.data['polar_path'].str.replace('/','\\')

        # Close file
        f.close()

        # Add airfoil objects
        self.airfoil_dir = os.path.join(os.path.dirname(path), "Airfoils")

        if os.path.isdir(self.airfoil_dir):
            for f in os.listdir(self.airfoil_dir):
                self.airfoils.append(Airfoil(os.path.join(self.airfoil_dir, f)))
        else:
            print(f"No airfoil folder found at {self.airfoil_dir}!")
            sys.exit(1)

        # Add polar objects
        self.polar_dir = os.path.join(os.path.dirname(path), "Polars")

        if os.path.isdir(self.polar_dir):
            for f in os.listdir(self.polar_dir):
                self.polars.append(Polar(os.path.join(self.polar_dir, f)))
        else:
            print(f"No polar folder found at {self.polar_dir}!")
            sys.exit(1)

if __name__ == "__main__":

    # Parse turbine file
    blade = Aero("data\\turbines\\DTU_10MW\\Aero\\DTU_10MW.bld")

    # Print turbine data
    print("Blade Name:", blade.name)
    print("Rotor Type:", blade.type)
    print("Inverted Foils:", blade.inverted_foils)
    print("Number of Blades:", blade.n_blades)
    print("Blade Data:")
    print(blade.data)
