""" 
Author:   T. Moreira da Fonte Fonseca Teles
Email:    tmoreiradafont@tudelft.nl
Date:     2025-02-07
License:  GNU GPL 3.0

Run the XFoil executable.

Classes:
    XFoil

Functions:
    None

Exceptions:
    None
"""

import os
import subprocess as sp
import sys


class XFoil:
    """
    A class to run the XFoil executable.
    
    Methods:
        __init__ -- spawn an XFoil process
        load_airfoil -- load an airfoil from a file
        run -- run XFoil at a given Re, M, and Alpha
    
    Attributes:
        path : str -- the path to the XFoil executable
        process : subprocess.Popen -- the XFoil process
    """

    def __init__(self, path):
        """
        Spawn an XFoil process.

        Parameters:
            path : str -- the path to the XFoil executable

        Returns:
            None
        """

        self.path = path

        if os.path.isfile(path):
            self.process = sp.Popen(path, stdin=sp.PIPE, stdout=sp.PIPE, stderr=sp.PIPE, text=True)
        else:
            print(f"No executable found at {path}!")
            sys.exit(1)

    def load_airfoil(self, path):
        """
        Load an airfoil from a file.

        Parameters:
            path : str -- the path to the airfoil file

        Returns:
            None
        """

        if os.path.isfile(path):
            self.process.stdin.write(f"LOAD {path}\n")
        else:
            print(f"No airfoil file found at {path}!")
            sys.exit(1)

    def run(self, re, mach, alpha, path, it=100):
        """
        Run XFoil at a given Re, M, and Alpha.

        Parameters:
            re : float -- the Reynolds number
            mach : float -- the Mach number
            alpha : float -- the angle of attack
            path : str -- the path to the output file
            it : int -- the maximum number of iterations (default 100)

        Returns:
            stdout -- the standard output stream
            stderr -- the standard error stream
        """

        # Create the output directory if it doesn't exist
        if not os.path.exists(os.path.dirname(path)):
            os.makedirs(os.path.dirname(path))

        # Queue XFoil commands
        self.process.stdin.write("OPER\n")
        self.process.stdin.write(f"ITER {it}\n")
        self.process.stdin.write(f"Visc {re}\n")
        self.process.stdin.write(f"Mach {mach}\n")
        self.process.stdin.write(f"ALFA {alpha}\n")
        self.process.stdin.write(f"DUMP {path}\n")
        self.process.stdin.write("\n")
        self.process.stdin.write("QUIT\n")

        # Run XFoil
        stdout, stderr = self.process.communicate()

        # Check for failed convergence
        if "Type \"!\" to continue iterating" in stdout:
            print("XFoil convergence failed")
            os.remove(path)
            sys.exit(1)

        return stdout, stderr
