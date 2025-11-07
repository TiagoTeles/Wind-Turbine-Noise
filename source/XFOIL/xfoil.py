""" 
Author:   T. Moreira da Fonte Fonseca Teles
Email:    tmoreiradafont@tudelft.nl
Date:     2025-11-07
License:  GNU GPL 3.0

Run the XFOIL executable.

Classes:
    XFoil

Functions:
    None

Exceptions:
    None
"""

import subprocess as sp
import sys

import numpy as np


class XFoil:
    """
    A class to run the XFOIL executable.

    Methods:
        __init__ -- initialise the XFoil class
        interpolate -- interpolate two airfoils
        run -- run XFOIL at a given Re, M, and Alpha

    Attributes:
        path : str -- path to the XFOIL executable
        cwd : str -- current working directory
        process : subprocess.Popen -- XFOIL process
    """

    def __init__(self, path, cwd):
        """
        Initialise the XFoil class.

        Parameters:
            path : str -- path to the XFOIL executable
            cwd : str -- current working directory

        Returns:
            None
        """

        self.path = path
        self.cwd = cwd

        # Create the XFOIL process
        self.process = sp.Popen(path, stdin=sp.PIPE, stdout=sp.PIPE, stderr=sp.PIPE, cwd=cwd,
                                text=True)

    def interpolate(self, path_0, path_1, path_out, fraction):
        """
        Interpolate two airfoils.

        Parameters:
            path_0 : str -- input airfoil file path
            path_1 : str -- input airfoil file path
            path_out : str -- output airfoil file path
            fraction : float -- interpolation fraction, [-]
        """

        # Check if the file paths are too long
        if len(path_0) > 64 or len(path_1) > 64 or len(path_out) > 64:
            print("File path is too long!")
            sys.exit(1)

        # Select the INTE environment
        self.process.stdin.write("INTE\n")

        # Load the first airfoil file
        self.process.stdin.write("F\n")
        self.process.stdin.write(f"{path_0}\n")

        # Load the second airfoil file
        self.process.stdin.write("F\n")
        self.process.stdin.write(f"{path_1}\n")

        # Set the interpolation fraction
        self.process.stdin.write(f"{fraction}\n")

        # Save the interpolated airfoil file
        self.process.stdin.write(f"interpolated_airfoil\n")
        self.process.stdin.write("PCOP\n")
        self.process.stdin.write(f"SAVE {path_out}\n")

        # Run the XFOIL commands
        stdout, _ = self.process.communicate()

        # Check for errors
        if "File OPEN error" in stdout:
            print("XFOIL file OPEN error!")
            sys.exit(1)

        elif "File READ error" in stdout:
            print("XFOIL file READ error!")
            sys.exit(1)

        elif "Output file exists" in stdout:
            print("Airfoil already exists! Skipping!")

    def run(self, path_in, path_out, re, mach, alpha, x_c_upper, x_c_lower, n_crit, max_iter):
        """
        Run XFOIL at a given Re, M, and Alpha.

        Parameters:
            path_in : str -- airfoil file path
            path_out : str -- dump file path
            re : float -- Reynolds number, [-]
            mach : float -- Mach number, [-]
            alpha : float -- angle of attack, [rad]
            x_c_upper : float -- transition point on the upper surface, [-]
            x_c_lower : float -- transition point on the lower surface, [-]
            n_crit : float -- critical amplification factor, [-]
            max_iter : int -- maximum number of XFOIL iterations, [-]
        """

        # Check if the file path is too long
        if len(path_in) > 64 or len(path_out) > 64:
            print("File path is too long!")
            sys.exit(1)

        # Convert alpha from [rad] to [deg]
        alpha = np.degrees(alpha)

        # Load the airfoil file
        self.process.stdin.write(f"LOAD {path_in}\n")

        # Set re, mach, and max_iter in the OPER environment
        self.process.stdin.write("OPER\n")
        self.process.stdin.write(f"Visc {re}\n")
        self.process.stdin.write(f"Mach {mach}\n")
        self.process.stdin.write(f"ITER {max_iter}\n")

        # Set x_c_upper, x_c_lower, and n_crit in the VPAR environment
        self.process.stdin.write("VPAR\n")
        self.process.stdin.write(f"XTR {x_c_upper} {x_c_lower}\n")
        self.process.stdin.write(f"N {n_crit}\n")
        self.process.stdin.write("\n")

        # Run the analysis at a given alpha
        self.process.stdin.write(f"Alfa {alpha}\n")

        # Save the boundary layer data
        self.process.stdin.write(f"DUMP {path_out}\n")
        self.process.stdin.write("\n")

        # Run the XFOIL commands
        stdout, _ = self.process.communicate()

        # Check for errors
        if "File OPEN error" in stdout:
            print("XFOIL file OPEN error!")
            sys.exit(1)

        elif "File READ error" in stdout:
            print("XFOIL file READ error!")
            sys.exit(1)

        elif "Type \"!\" to continue iterating" in stdout:
            print("XFOIL convergence failed!")
            sys.exit(1)
