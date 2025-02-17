import os
import sys
import shutil

import numpy as np
import pandas as pd

from qblade import QBlade
from turbine import Turbine

# QBlade configuration
DLL_PATH = "bin\\QBlade\\QBlade.dll"    # DLL path
CL_DEVICE = 1                           # OpenCL device
GROUP_SIZE = 32                         # OpenCL work-group size

# Simulation specification
TURBINE_PATH = "data\\turbines\\DTU_10MW\\DTU_10MW_RWT.trb" # Turbine path
SIMULATION_PATH = "data\\simulations\\test.sim"             # Turbine name
N_TIMESTEP = 500                                            # Number of timesteps


if not os.path.exists("temp"):

    os.makedirs("temp")

    shutil.copyfile(SIMULATION_PATH, "temp\\simulation.sim")
    shutil.copytree("data\\turbines\\DTU_10MW", "temp\\DTU_10MW_RWT")

    # Ensure a matching file is found
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
    qblade.loadSimDefinition("temp\\simulation.sim".encode("utf-8"))

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

results = pd.read_csv("temp\\results.txt", skiprows=2, delimiter="\t").tail(1)

import matplotlib.pyplot as plt


for i in range(20, 29, 2):
    print(i)
    r = results[f"Radius_Blade_1_PAN_{i}_[m]"]
    U = results[f"Total_Velocity_Blade_1_PAN_{i}_[m/s]"]

    M = np.array(U) / 340
    H = 100
    L = min(0.7*H, 42)
    b =  results[f"Radius_Blade_1_PAN_{i+1}_[m]"] - results[f"Radius_Blade_1_PAN_{i}_[m]"]
    z = 100
    I = 0.01
    K_e = 0.75/L

    spl = []


    for f in np.logspace(0, 4, 50):

        omega = 2 * np.pi * f
        K_x = omega / U

        K_x_hat = K_x / K_e


        c = [5.38, 5.38, 5.38, 5.38, 5.3886, 5.4212, 5.4865, 5.5887, 5.7247, 5.8817, 6.0346, 6.1478, 6.202, 6.195, 6.1292, 6.0096, 5.8432, 5.64, 5.4107, 5.1613, 4.8974, 4.6255, 4.3519, 4.0827, 3.822, 3.5724, 3.3364, 3.1161, 2.913, 2.7275, 2.5595, 2.4087, 2.266, 2.1175, 1.9588, 1.7913, 1.6013, 1.3858, 1.1384, 0.83354]
        K_x_line = K_x * c[i] / 2
        beta = np.sqrt(1-M*M)
        S = 1 / (2*np.pi * K_x_line/np.square(beta) + 1/(1+2.4*K_x_line/np.square(beta)))

        SPL_H = 10 * np.log10(np.pow(M, 5) * L * b / (2*z*z) *np.square(I) * np.square(1.225) * np.pow(340, 4) * np.pow(K_x_hat, 3) / (np.pow((1+np.square(K_x_hat)), 7/3))) + 58.4
        LFC = 10 * np.square(S) * M  * np.square(K_x_line) / np.square(beta)
        SPL_amiet = SPL_H + 10 * np.log10(LFC/(1+LFC))

        spl.append(SPL_amiet)

    plt.plot(np.logspace(0, 4, 50), spl)
    plt.xscale('log')
    plt.grid()

plt.show()



