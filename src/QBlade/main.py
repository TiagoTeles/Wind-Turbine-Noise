import os
import sys
# from ctypes import *
from QBladeLibrary import QBladeLibrary
dll_name = 'bin/QBlade/'
turbine_path = 'data/QBlade/DTU_10MW/'
turbine_name = 'DTU_10MW'

# On Windows systems, we update the PATH environment variable to include the QBlade directory.
# This ensures that the required SSL libraries (e.g., libssl and libcrypto) are properly located and loaded.
# If experiencing issues with this DLL in a Windows Python environment see:
# https://docs.qblade.org/src/license/license_files.html#resolving-openssl-issues-on-windows
if os.name == 'nt':  # 'nt' indicates Windows
    os.environ["PATH"] = os.path.abspath(dll_name) + ";" + os.environ.get("PATH", "")

# Search the directory below for library files matching the pattern QBlade*.dll or QBlade*.so
dll_files = []

for f in os.listdir(dll_name):
    if 'QBlade' in f and ('.dll' in f or '.so' in f):
        dll_files.append(f)

# Check if any matching files are found
if not dll_files:
    print('No matching QBlade*.dll or QBlade*.so files found in the specified directory:', os.path.abspath(dll_name))
    sys.exit(1)  # Exit the script with a non-zero status to indicate an error

# Use the first matching file
dll_file_path = os.path.join(dll_name, dll_files[0])

        except Exception as e:
            raise RuntimeError(f"Could not load the library at {self.lib_path}: {e}") from e

# Create an object of the class 'QBladeLibrary' that contains the API
QBlade = QBladeLibrary(dll_file_path)

# Creation of a QBlade instance from the library
QBlade.createInstance(1, 32)

# Loading a project or sim-file, in this case the DTU_10MW_Demo project or simulation definition file
QBlade.loadSimDefinition(os.path.join(turbine_path, turbine_name, '.sim'))

# Initializing the sim and ramp-up phase, call before starting the simulation loop
QBlade.initializeSimulation()

    def unload_library(self):
        """Close the QBlade instance if it exists."""

        # Close QBlade instance
        try:
            self.closeInstance()
            print("QBlade instance closed.")

    #advance the simulation
    success = QBlade.advanceTurbineSimulation()

    # Check if the simulation step was successful
    if not success:  # If success is False, exit the loop
        print(f"Simulation failed at timestep {i}. Exiting loop.")
        break

# Storing the simulation results in QBlade ASCII format in the file NREL_5MW_Sample_results.txt
QBlade.exportResults(0, turbine_path, turbine_name, b"")

# Unloading the qblade library
QBlade.unload_library() 
