""" 
Author:   T. Moreira da Fonte Fonseca Teles
Email:    tmoreiradafont@tudelft.nl
Date:     2025-02-14
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
        name : str -- name of the airfoil object
        coordinates : pd.DataFrame -- coordinates of the airfoil [-]
        path : str -- path to the .afl file
    """

    def __init__(self, path):
        """
        Parse the airfoil file.

        Arguments:
            path : str -- path to the .afl file

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
        self.name = f.readline().strip("\n")

        # Read x and y coordinates
        self.coordinates = pd.read_csv(f, names=['x', 'y'], sep=r"\s+")

        # Close file
        f.close()

if __name__ == "__main__":

    # Parse file
    airfoil = Airfoil("data\\turbines\\DTU_10MW\\Aero\\Airfoils\\FFA_W3_241.afl")

    # Print contents
    print("Airfoil Name:", airfoil.name)
    print("Airfoil X Coordinates: [-]")
    print(airfoil.coordinates['x'])
    print("Airfoil Y Coordinates: [-]")
    print(airfoil.coordinates['y'])
