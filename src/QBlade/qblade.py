import os
import sys
from ctypes import *
from QBladeLibrary import QBladeLibrary

# Define the directory where the QBlade library is located
dll_directory = "bin/QBlade/"

# On Windows systems, we update the PATH environment variable to include the QBlade directory. 
# This ensures that the required SSL libraries (e.g., libssl and libcrypto) are properly located and loaded.
# If experiencing issues with this DLL in a Windows Python environment see:
# https://docs.qblade.org/src/license/license_files.html#resolving-openssl-issues-on-windows
if os.name == 'nt':  # 'nt' indicates Windows
    os.environ["PATH"] = os.path.abspath(dll_directory) + ";" + os.environ.get("PATH", "")

# Search the directory below for library files matching the pattern QBlade*.dll or QBlade*.so
dll_files = [f for f in os.listdir(dll_directory) if 'QBlade' in f and ('.dll' in f or '.so' in f)]

# Check if any matching files are found
if not dll_files:
    print('No matching QBlade*.dll or QBlade*.so files found in the specified directory:',os.path.abspath(dll_directory))
    sys.exit(1)  # Exit the script with a non-zero status to indicate an error

# Use the first matching file
dll_file_path = os.path.join(dll_directory, dll_files[0])

# Display the selected shared library file
print(f'Using shared library file: {dll_file_path}')

# Create an object of the class 'QBladeLibrary' that contains the API
QBLADE = QBladeLibrary(dll_file_path)    

# Creation of a QBlade instance from the library
QBLADE.createInstance(1,32)

# Loading a project or sim-file, in this case the DTU_10MW_Demo project or simulation definition file
#QBLADE.loadSimDefinition(b"./DTU_10MW_Demo.sim") #uncomment this line to load a simulation definition file
QBLADE.loadProject(b"data/QBlade/NREL_5MW.qpr") 

# Initializing the sim and ramp-up phase, call before starting the simulation loop
QBLADE.initializeSimulation()

# We will run the simulation for 500 steps before storing the results
number_of_timesteps = 500

# Start of the simulation loop
for i in range(number_of_timesteps):

    #advance the simulation
    success = QBLADE.advanceTurbineSimulation() 
    
    # Check if the simulation step was successful
    if not success:  # If success is False, exit the loop
        print(f"Simulation failed at timestep {i}. Exiting loop.")
        break

# Storing the finished simulation in a project as NREL_5MW_Sample_completed, you can open this file to view the results of the simulation inside QBlade's GUI
QBLADE.storeProject(b"./NREL_5MW_Sample_completed.qpr")

# Storing the simulation results in QBlade ASCII format in the file NREL_5MW_Sample_results.txt
QBLADE.exportResults(0,b"./",b"NREL_5MW_Sample_results",b"")

# Unloading the qblade library
QBLADE.unload() 