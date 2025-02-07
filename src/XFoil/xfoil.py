import os
import subprocess as sp
import sys


class XFoil:
    """A class to interact with XFoil."""

    def __init__(self, path):
        """
        Spawn an XFoil process.

        Keyword arguments:
            path -- the path to the XFoil executable
        """

        if os.path.isfile(path):
            self.process = sp.Popen(path, stdin=sp.PIPE, stdout=sp.PIPE, stderr=sp.PIPE, text=True)
        else:
            print(f"No executable file found at {path}")
            sys.exit(1)

    def load_airfoil(self, path):
        """
        Load an airfoil from a file.

        Keyword arguments:
            path -- the path to the airfoil file
        """

        if os.path.isfile(path):
            self.process.stdin.write(f"LOAD {path}\n")
        else:
            print(f"No airfoil file found at {path}")
            sys.exit(1)

    def run(self, re, mach, alpha, path):
        """
        Run XFoil at a given Re, M, and Alpha.
        
        Keyword arguments:
            re -- the Reynolds number
            mach -- the Mach number
            alpha -- the angle of attack
            path -- the path to the output file

        Returns:
            stdout -- the standard output stream
            stderr -- the standard error stream
        """

        # Create the output directory if it doesn't exist
        if not os.path.exists(os.path.dirname(path)):
            os.makedirs(os.path.dirname(path))

        # Queue XFoil commands
        self.process.stdin.write("OPER\n")
        self.process.stdin.write(f"Visc {re}\n")
        self.process.stdin.write(f"Mach {mach}\n")
        self.process.stdin.write(f"ALFA {alpha}\n")
        self.process.stdin.write(f"DUMP {path}\n")

        return self.process.communicate()
