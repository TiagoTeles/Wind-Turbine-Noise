import os
import sys
import shutil

from qblade import QBlade
from simulation import Simulation

# QBlade configuration
DLL_PATH = "bin\\QBlade\\QBlade.dll"    # DLL path
CL_DEVICE = 1                           # OpenCL device
GROUP_SIZE = 32                         # OpenCL work-group size

# Simulation specification
TURBINE_NAME = "DTU_10MW"                                   # Turbine name
TURBINE_PATH = "data\\turbines\\DTU_10MW\\DTU_10MW_RWT.trb" # Turbine definition path
SIMULATION_PATH = "data\\simulations\\DTU_10MW_RWT.sim"     # Simulation definition path
WORKING_PATH = "temp"                                       # Working directory path
N_TIMESTEP = 500                                            # Number of timesteps

# Copy files to working directory
shutil.copytree(os.path.dirname(TURBINE_PATH), os.path.join(WORKING_PATH, TURBINE_NAME))
shutil.copyfile(SIMULATION_PATH, os.path.join(WORKING_PATH, TURBINE_NAME + ".sim"))

# Read files
simulation = Simulation(os.path.join(WORKING_PATH, TURBINE_NAME + ".sim"))

# Start QBlade
if os.path.exists(DLL_PATH):
    print(f"Using shared library file: {DLL_PATH}!")
else:
    print(f"No shared library file was found in {DLL_PATH}!")
    sys.exit(1)

# Create an object that contains the API
qblade = QBlade(DLL_PATH)

# Create a QBlade instance from the library
qblade.createInstance(CL_DEVICE, GROUP_SIZE)

# Load a simulation definition file
qblade.loadSimDefinition(os.path.join(WORKING_PATH, TURBINE_NAME + ".sim").encode("utf-8"))

# Initializing the simulation
qblade.initializeSimulation()

# Run the simulation
for i in range(N_TIMESTEP):

    # Advance the simulation one timestep
    success = qblade.advanceTurbineSimulation()

    # Ensure the simulation step was successful
    if not success:
        print(f"Simulation failed at timestep {i}!")
        break

# Store the simulation results
qblade.exportResults(0, "temp".encode("utf-8"), "results".encode("utf-8"), b"")

# Unload the QBlade library
qblade.unload_library()
