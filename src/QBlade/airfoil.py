""" 
Author:   T. Moreira da Fonte Fonseca Teles
Email:    tmoreiradafont@tudelft.nl
Date:     2025-02-07
License:  GNU GPL 3.0

Store airfoil data.

Classes:
    Airfoil

Functions:
    None

Exceptions:
    None
"""

import os
import sys
import pandas as pd


class Airfoil:
    """ 
    
    A class to store the airfoil data. 
    
    Methods:
        __init__ -- parse the airfoil file
        
    Attributes:
        path : str -- the path to the .afl file
        airfoil_name : str -- the name of the airfoil
        coordinates : pd.DataFrame -- the x and y coordinates of the airfoil
    """

    def __init__(self, path):
        """
        Parse the airfoil file.

        Arguments:
            path : str -- the path to the .afl file

        Returns:
            None
        """

        self.path = path

        # Open file
        if os.path.isfile(path):
            f = open(path, "r", encoding="utf-8")
        else:
            print(f"No airfoil file found at {path}!")
            sys.exit(1)

        # Parse data in file
        self.airfoil_name = f.readline().strip("\n")

        # Read x and y coordinates
        self.coordinates = pd.read_csv(f, names=['x', 'y'], sep=r"\s+")

        # Close file
        f.close()
