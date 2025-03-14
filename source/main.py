"""
Author:   T. Moreira da Fonte Fonseca Teles
Email:    tmoreiradafont@tudelft.nl
Date:     2025-03-14
License:  GNU GPL 3.0

Main script to run a QBlade simulation using the Python API.

Classes:
    None

Functions:
    None

Exceptions:
    None
"""

import os
import sys
import shutil

from QBlade.qblade import QBlade
from QBlade.simulation import Simulation

# QBlade configuration
QBLADE_PATH = "bin\\QBlade\\QBladeCE_2.0.8.5.dll"                   # DLL path
QBLADE_CL_DEVICE = 1                                                # OpenCL device
QBLADE_GROUP_SIZE = 32                                              # OpenCL work-group size

# Simulation specification
TURBINE_PATH = "data\\turbines\\DTU_10MW_v1.3\\DTU_10MW_RWT.trb"    # Turbine definition path
SIMULATION_PATH = "data\\simulations\\DTU_10MW_RWT.sim"             # Simulation definition path
WORKING_DIR = "temp"                                                # Working directory
N_TIMESTEP = 500                                                    # Number of timesteps, [-]

# Copy files to the working directory
turbine_dir = os.path.dirname(TURBINE_PATH)
turbine_name = os.path.splitext(os.path.basename(TURBINE_PATH))[0]

simulation_name = os.path.splitext(os.path.basename(SIMULATION_PATH))[0]

turbine_work_dir = os.path.join(WORKING_DIR, simulation_name, turbine_name)
simulation_work_path = os.path.join(WORKING_DIR, simulation_name, simulation_name + ".sim")

if not os.path.exists(WORKING_DIR):
    shutil.copytree(turbine_dir, turbine_work_dir)
    shutil.copyfile(SIMULATION_PATH, simulation_work_path)

# Start QBlade
if os.path.exists(QBLADE_PATH):
    print(f"Using shared library file: {QBLADE_PATH}!")
else:
    print(f"No shared library file was found in {QBLADE_PATH}!")
    sys.exit(1)

qblade = QBlade(QBLADE_PATH)
qblade.createInstance(QBLADE_CL_DEVICE, QBLADE_GROUP_SIZE)

# Configure the simulation
simulation = Simulation(simulation_work_path)

simulation.turbine.write("BEMSPEEDUP", 100.0)   # Speedup convergence
simulation.turbine.write("STRUCTURALFILE", "")  # Disable structural simulation
simulation.turbine.write("CONTROLLERFILE", "")  # Disable controller simulation
simulation.turbine.write("PARAMETERFILE", "")   # Disable controller simulation
simulation.turbine.write("CONTROLLERTYPE", 0)   # Disable controller simulation

# Setup the simulation
qblade.loadSimDefinition(simulation_work_path.encode("utf-8"))
qblade.initializeSimulation()

# Run the simulation
for i in range(N_TIMESTEP):

    # Advance the simulation one timestep
    success = qblade.advanceTurbineSimulation()

    # Ensure the simulation step was successful
    if not success:
        print(f"Simulation failed at timestep {i}!")
        sys.exit(1)

# Save the simulation results
qblade.exportResults(0, os.path.join(WORKING_DIR, simulation_name).encode("utf-8"), simulation_name.encode("utf-8"), b"")

# Unload the QBlade library
qblade.unload_library()
