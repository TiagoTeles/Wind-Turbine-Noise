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

import numpy as np
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
        path : str -- path to the XFoil executable
        process : subprocess.Popen -- XFoil process
    """

    def __init__(self, path):
        """
        Spawn an XFoil process.

        Arguments:
            path : str -- path to the XFoil executable

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

        Arguments:
            path : str -- path to the airfoil file

        Returns:
            None
        """

        if os.path.isfile(path):
            self.process.stdin.write(f"LOAD {path}\n")
        else:
            print(f"No airfoil file found at {path}!")
            sys.exit(1)

    def run(self, re, mach, alpha, path, xtr_top=1.0, xtr_bot=1.0, n_crit=9.0, it=10):
        """
        Run XFoil at a given Re, M, and Alpha.

        Arguments:
            re : float -- Reynolds number, [-]
            mach : float -- Mach number, [-]
            alpha : float -- angle of attack, [rad]
            path : str -- path to the output file
            xtr_top : float -- transition point on the top surface (default 1.0), [-]
            xtr_bot : float -- transition point on the bottom surface (default 1.0), [-]
            n_crit : float -- critical amplification factor (default 9.0), [-]
            it : int -- maximum number of iterations (default 100), [-]

        Returns:
            stdout -- the standard output stream
            stderr -- the standard error stream
        """

        # Set Re, Ma, and ITER in OPER environment
        self.process.stdin.write("OPER\n")
        self.process.stdin.write(f"Visc {re}\n")
        self.process.stdin.write(f"Mach {mach}\n")
        self.process.stdin.write(f"ITER {it}\n")

        # Set Xtr_top, Xtr_bot, and N_crit in VPAR environment
        self.process.stdin.write("VPAR\n")
        self.process.stdin.write(f"XTR {xtr_top} {xtr_bot}\n")
        self.process.stdin.write(f"N {n_crit}\n")
        self.process.stdin.write("\n")

        # Run analysis at a given alpha
        self.process.stdin.write(f"Alfa {np.radians(alpha)}\n")

        # Save results to a file
        if not os.path.exists(os.path.dirname(path)):
            os.makedirs(os.path.dirname(path))

        self.process.stdin.write(f"DUMP {path}\n")
        self.process.stdin.write("\n")

        # Quit XFoil
        self.process.stdin.write("QUIT\n")

        # Run XFoil commands
        stdout, stderr = self.process.communicate()

        # Check for failed convergence
        if "Type \"!\" to continue iterating" in stdout:
            print("XFoil convergence failed!")
            os.remove(path)
            sys.exit(1)

        return stdout, stderr
