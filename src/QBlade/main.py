import os
import QBladeLib
import sys

# QBlade configuration
DLL_PATH = "bin/QBlade/"    # DLL path
DLL_NAME = "QBlade.dll"     # DLL name
CL_DEVICE = 1               # OpenCL device
GROUP_SIZE = 32             # OpenCL work-group size

# Simulation specification
TURBINE_PATH = "data/QBlade/DTU_10MW/"  # Turbine path
TURBINE_NAME = "DTU_10MW.sim"           # Turbine name
N_TIMESTEP = 500                        # Number of timesteps

# Search DLL_PATH for files matching DLL_NAME
dll_files = []

for f in os.listdir(DLL_PATH):
    if f == DLL_NAME:
        dll_files.append(f)

# Ensure a matching file is found
if dll_files:
    print(f"Using shared library file: {os.path.join(DLL_PATH, DLL_NAME)}!")
else:
    print(f"No shared library file named {DLL_NAME} was found in {DLL_PATH}!")
    sys.exit(1)

# Create an object that contains the API
QBlade = QBladeLib.QBlade(os.path.join(DLL_PATH, DLL_NAME))

# Create a QBlade instance from the library
QBlade.createInstance(CL_DEVICE, GROUP_SIZE)

# Load a simulation definition file
QBlade.loadSimDefinition(os.path.join(TURBINE_PATH, TURBINE_NAME).encode("utf-8"))

# Initializing the simulation
QBlade.initializeSimulation()

# Run the simulation
for i in range(N_TIMESTEP):

    # Advance the simulation one timestep
    success = QBlade.advanceTurbineSimulation()

    # Ensure the simulation step was successful
    if not success:
        print(f"Simulation failed at timestep {i}!")
        break

# Store the simulation results
QBlade.exportResults(0, TURBINE_PATH.encode("utf-8"), os.path.splitext(TURBINE_NAME)[0].encode("utf-8"), b"")

# Unload the QBlade library
QBlade.unload_library()
