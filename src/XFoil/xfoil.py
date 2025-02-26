""" 
Author:   T. Moreira da Fonte Fonseca Teles
Email:    tmoreiradafont@tudelft.nl
Date:     2025-02-26
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

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd


class XFoil:
    """
    A class to run the XFoil executable.

    Methods:
        __init__ -- spawn an XFoil process
        run -- run XFoil at a given Re, M, and Alpha

    Attributes:
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

        if os.path.isfile(path):
            self.process = sp.Popen(path, stdin=sp.PIPE, stdout=sp.PIPE, stderr=sp.PIPE, text=True)
        else:
            print(f"No executable found at {path}!")
            sys.exit(1)

    def run(self, path, re, mach, alpha, xtr_top=1.0, xtr_bot=1.0, n_crit=9.0, it=10):
        """
        Run XFoil at a given Re, M, and Alpha.

        Arguments:
            path : str -- airfoil file path
            re : float -- Reynolds number, [-]
            mach : float -- Mach number, [-]
            alpha : float -- angle of attack, [rad]
            xtr_top : float -- transition point on the top surface (default 1.0), [-]
            xtr_bot : float -- transition point on the bottom surface (default 1.0), [-]
            n_crit : float -- critical amplification factor (default 9.0), [-]
            it : int -- maximum number of iterations (default 10), [-]

        Returns:
            bl_top : pandas.DataFrame -- boundary layer data on the top surface
            bl_bot : pandas.DataFrame -- boundary layer data on the bottom surface
        """

        # Load airfoil file
        if os.path.isfile(path):
            self.process.stdin.write(f"LOAD {path}\n")
        else:
            print(f"No airfoil file found at {path}!")
            sys.exit(1)

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
        self.process.stdin.write(f"Alfa {np.degrees(alpha)}\n")

        # Save results to a file
        output_path = os.path.join("tmp\\XFoil", os.path.basename(path))

        if not os.path.exists(os.path.dirname(output_path)):
            os.makedirs(os.path.dirname(output_path))

        self.process.stdin.write(f"DUMP {output_path}\n")
        self.process.stdin.write("\n")

        # Run XFoil commands
        stdout, _ = self.process.communicate()

        # Check for failed convergence
        if "Type \"!\" to continue iterating" in stdout:
            print("XFoil convergence failed!")
            print(f"Airfoil: {path}, Re: {re}, Ma: {mach}, Alpha: {alpha}")
            sys.exit(1)

        # Read output file
        data = pd.read_csv(output_path, sep=r"\s+", skiprows=1,
                           names=["s/c", "x/c", "y/c", "U_e/U", "delta_star", "theta", "C_f", "H"])

        # Remove output file
        os.remove(output_path)

        # Filter and sort data
        le_index = data.idxmin()["x/c"]
        bl_top = data[data["x/c"] <= 1.0].iloc[:le_index+1][::-1]
        bl_bot = data[data["x/c"] <= 1.0].iloc[le_index:]

        return bl_top, bl_bot

# if __name__ == "__main__":

#     XFOIL_PATH = "bin\\XFoil\\xfoil.exe"
#     AIRFOIL_PATH = "data\\turbines\\DTU_10MW\\Aero\\Airfoils\\FFA_W3_241.afl"

#     xfoil = XFoil(XFOIL_PATH)
#     top, bot = xfoil.run(AIRFOIL_PATH, 1E6, 0.2, 0)

#     plt.plot(top["x/c"], top["y/c"])
#     plt.plot(bot["x/c"], bot["y/c"])
#     plt.show()

#     plt.plot(top["x/c"], top["U_e/U"])
#     plt.plot(bot["x/c"], bot["U_e/U"])
#     plt.show()

#     plt.plot(top["x/c"], top["delta_star"])
#     plt.plot(bot["x/c"], bot["delta_star"])
#     plt.show()

#     plt.plot(top["x/c"], top["theta"])
#     plt.plot(bot["x/c"], bot["theta"])
#     plt.show()

#     plt.plot(top["x/c"], top["C_f"])
#     plt.plot(bot["x/c"], bot["C_f"])
#     plt.show()

#     plt.plot(top["x/c"], top["H"])
#     plt.plot(bot["x/c"], bot["H"])
#     plt.show()
