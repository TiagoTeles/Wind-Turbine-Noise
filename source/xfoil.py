""" 
Author:   T. Moreira da Fonte Fonseca Teles
Email:    tmoreiradafont@tudelft.nl
Date:     2025-03-17
License:  GNU GPL 3.0

Run the XFOIL executable.

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

import numpy as np
import pandas as pd


class XFoil:
    """
    A class to run the XFOIL executable.

    Methods:
        __init__ -- spawn an XFoil process
        interpolate -- interpolate between two airfoil files
        run -- run XFOIL at a given Re, M, and Alpha

    Attributes:
        cwd : str -- working directory
        path : str -- path to the XFOIL executable
        process : subprocess.Popen -- XFOIL process
    """

    def __init__(self, path, cwd):
        """
        Spawn an XFoil process.

        Arguments:
            path : str -- path to the XFOIL executable
            cwd : str -- working directory

        Returns:
            None
        """

        # Initialise the attributes
        self.path = path
        self.cwd = cwd

        # Create the XFOIL process
        self.process = sp.Popen(path, stdin=sp.PIPE, stdout=sp.PIPE, stderr=sp.PIPE, cwd=cwd, 
                                text=True)

    def interpolate(self, path_1, path_2, fraction):
        """
        Interpolate between two airfoil files.

        Arguments:
            path_1 : str -- input airfoil 1 file path
            path_2 : str -- input airfoil 2 file path
            fraction : float -- interpolation fraction, [-]

        Returns:
            path_out : str -- interpolated airfoil file path
        """

        # Select the INTE environment
        self.process.stdin.write("INTE\n")

        # Load the first airfoil file
        self.process.stdin.write("F\n")
        self.process.stdin.write(f"{path_1}\n")

        # Load the second airfoil file
        self.process.stdin.write("F\n")
        self.process.stdin.write(f"{path_2}\n")

        # Set the interpolation fraction
        self.process.stdin.write(f"{fraction}\n")

        # Determine the output file name and path
        name_1 = os.path.splitext(os.path.basename(path_1))[0]
        name_2 = os.path.splitext(os.path.basename(path_2))[0]

        name_out = f"{name_1}_f{fraction:.2f}"
        path_out = os.path.join("Interpolated", name_out + ".afl")

        # Save the interpolated airfoil file
        if not os.path.exists(os.path.join(self.cwd, "Interpolated")):
            os.makedirs(os.path.join(self.cwd, "Interpolated"))

        self.process.stdin.write(f"{name_out}\n")
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

        elif "Output file exists." in stdout:
            print("Airfoil already exists! Skipping!")

        return os.path.join(self.cwd, path_out)

    def run(self, path, re, mach, alpha, xtr_top=1.0, xtr_bot=1.0, n_crit=9.0, it=10):
        """
        Run XFOIL at a given Re, M, and Alpha.

        Arguments:
            path : str -- airfoil file path
            re : float -- Reynolds number, [-]
            mach : float -- Mach number, [-]
            alpha : float -- angle of attack, [rad]
            xtr_top : float -- transition point on the top surface, [-]
            xtr_bot : float -- transition point on the bottom surface, [-]
            n_crit : float -- critical amplification factor, [-]
            it : int -- maximum number of iterations, [-]

        Returns:
            bl_top : pandas.DataFrame -- boundary layer data on the top surface
            bl_bot : pandas.DataFrame -- boundary layer data on the bottom surface
        """

        # Load the airfoil file
        self.process.stdin.write(f"LOAD {path}\n")

        # Set the Re, Ma, and ITER in the OPER environment
        self.process.stdin.write("OPER\n")
        self.process.stdin.write(f"Visc {re}\n")
        self.process.stdin.write(f"Mach {mach}\n")
        self.process.stdin.write(f"ITER {it}\n")

        # Set the Xtr_top, Xtr_bot, and N_crit in the VPAR environment
        self.process.stdin.write("VPAR\n")
        self.process.stdin.write(f"XTR {xtr_top} {xtr_bot}\n")
        self.process.stdin.write(f"N {n_crit}\n")
        self.process.stdin.write("\n")

        # Run the analysis at a given alpha
        self.process.stdin.write(f"Alfa {np.degrees(alpha)}\n")

        # Determine the output file name and path
        name = os.path.splitext(os.path.basename(path))[0]

        name_out = f"Re{re:.2E}_Ma{mach:.2f}_AOA{alpha:.2f}"
        path_out = os.path.join("Boundary_Layer", name_out + ".bl")

        # Save the results to a file
        if not os.path.exists(os.path.join(self.cwd, "Boundary_Layer")):
            os.makedirs(os.path.join(self.cwd, "Boundary_Layer"))

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

        # Read the output file
        data = pd.read_csv(os.path.join(self.cwd, path_out), sep=r"\s+", skiprows=1,
                           names=["s/c", "x/c", "y/c", "U_e/U", "delta_star", "theta", "C_f", "H"])

        # Filter and sort the boundary layer data
        le_index = data["x/c"].idxmin()
        bl_top = data[data["x/c"] <= 1.0].iloc[:le_index + 1][::-1]
        bl_bot = data[data["x/c"] <= 1.0].iloc[le_index:]

        return bl_top, bl_bot
