""" 
Author:   T. Moreira da Fonte Fonseca Teles
Email:    tmoreiradafont@tudelft.nl
Date:     2025-02-07
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
        path : str -- the path to the .plr file
        name : str -- the name of the polar
        airfoil_path : str -- the path to the .afl file
        airfoil_thickness : float -- the thickness of the airfoil
        is_decomposed : bool -- whether the polar is decomposed
        reynolds : float -- the Reynolds number
        data : pd.DataFrame -- the AOA, CL, CD, and CM data
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
            print(f"No polar file found at {path}!")
            sys.exit(1)

        # Parse data in the file
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
