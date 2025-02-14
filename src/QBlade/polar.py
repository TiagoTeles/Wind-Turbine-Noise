""" 
Author:   T. Moreira da Fonte Fonseca Teles
Email:    tmoreiradafont@tudelft.nl
Date:     2025-02-14
License:  GNU GPL 3.0

Store polar data.

Classes:
    Polar

Functions:
    None

Exceptions:
    None
"""

import os
import sys
import pandas as pd

from misc import parse

class Polar:
    """ 
    A class to store the polar data.

    Methods:
        __init__ -- parse the polar file

    Attributes:
        airfoil_path : str -- the path to the .afl file
        airfoil_thickness : float -- the thickness of the airfoil
        data : pd.DataFrame -- the AOA, CL, CD, and CM data
        is_decomposed : bool -- whether the polar is decomposed
        name : str -- the name of the polar
        path : str -- the path to the .plr file
        reynolds : float -- the Reynolds number
    """

    def __init__(self, path):
        """
        Parse the polar file.

        Arguments:
            path : str -- the path to the .plr file

        Returns:
            None
        """

        self.path = path

        # Open file
        if os.path.isfile(path):
            f = open(path, "r", encoding="utf-8")
        else:
            print(f"No file found at {path}!")
            sys.exit(1)

        # Parse data in file
        self.name              = parse(f,    "POLARNAME", 0,   str)
        self.airfoil_path      = parse(f,     "FOILNAME", 0,   str)
        self.airfoil_thickness = parse(f,    "THICKNESS", 0, float)
        self.is_decomposed     = parse(f, "ISDECOMPOSED", 0,  bool)
        self.reynolds          = parse(f,     "REYNOLDS", 1, float)

        # Format parsed data
        self.airfoil_path = os.path.join(os.path.dirname(path), self.airfoil_path.replace("/", "\\"))
        self.airfoil_thickness /= 100.0

        # Read AOA, CL, CD, and CM
        self.data = pd.read_csv(f, delimiter=r"\s+", skiprows=2)

        # Close file
        f.close()

if __name__ == "__main__":

    # Parse File
    polar = Polar("data\\turbines\\DTU_10MW\\Aero\\Polars\\FFA_W3_241_t24.1_dtu_10mw_Polar_RE1.00E+06.plr")

    # Print Contents
    print("Polar Name:", polar.name)
    print("Airfoil Path:", polar.airfoil_path)
    print("Airfoil Thickness:", polar.airfoil_thickness)
    print("is Decomposed:", polar.is_decomposed)
    print("Reynolds Number:", polar.reynolds)
    print("Polar AOA:")
    print(polar.data['AOA'])
    print("Polar CL:")
    print(polar.data['CL'])
    print("Polar CD:")
    print(polar.data['CD'])
    print("Polar CM:")
    print(polar.data['CM'])
