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

    def run(self, re, mach, alpha, path, max_iter=1000):
        """
        Run XFoil at a given Re, M, and Alpha.

        Keyword arguments:
            re -- the Reynolds number
            mach -- the Mach number
            alpha -- the angle of attack
            path -- the path to the output file
            max_iter -- the maximum number of iterations (default 100)

        Returns:
            stdout -- the standard output stream
            stderr -- the standard error stream
        """

        # Create the output directory if it doesn't exist
        if not os.path.exists(os.path.dirname(path)):
            os.makedirs(os.path.dirname(path))

        # Queue XFoil commands
        self.process.stdin.write("OPER\n")
        self.process.stdin.write(f"ITER {max_iter}\n")
        self.process.stdin.write(f"Visc {re}\n")
        self.process.stdin.write(f"Mach {mach}\n")
        self.process.stdin.write(f"ALFA {alpha}\n")
        self.process.stdin.write(f"DUMP {path}\n")
        self.process.stdin.write("\n")
        self.process.stdin.write("QUIT\n")

        # Run XFoil commands
        stdout, stderr = self.process.communicate()

        # Check for convergence problems
        if "Type \"!\" to continue iterating" in stdout:
            print("VISCAL: Convergence failed")
            os.remove(path)
            # sys.exit(1)

        return stdout, stderr

for file in os.listdir("data/QBlade/DTU_10MW/Aero/Airfoils"):

    filepath = os.path.join("data/QBlade/DTU_10MW/Aero/Airfoils", file)

    test = XFoil("bin/XFoil/XFoil.exe")
    test.load_airfoil(filepath)
    test.run(1E6, 0.2, 0.0, filepath.replace(".afl", ".dat").replace("Airfoils", "Boundaries"))
