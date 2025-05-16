"""
Author:   T. Moreira da Fonte Fonseca Teles
Email:    tmoreiradafont@tudelft.nl
Date:     2025-05-14
License:  GNU GPL 3.0
"""

import os

import numpy as np

from aeroacoustics.trailing_edge_noise import te_noise
from coordinates import freestream_to_airfoil, airfoil_to_blade, blade_to_hub, hub_to_nacelle, nacelle_to_turbine
from miscellaneous import octave, surface_roughness_length, turbulence_intensity, turbulence_length_scale
from plot.aeroacoustic import plot_spectrum
from QBlade.simulation import Simulation
from settings import *


# Create the QBlade objects
print("Reading the QBlade objects...")

simulation = Simulation(SIMULATION_PATH)

# Run the simulation
print("Running the QBlade simulation...")

if not os.path.exists(SIMULATION_PATH.replace(".sim", ".txt")):

    # Initialise the simulation
    simulation.initialise(QBLADE_PATH, QBLADE_CL_DEVICE, QBLADE_GROUP_SIZE)

    # Run the simulation
    simulation.run(N_TIMESTEP)

    # Close the simulation
    simulation.close()

else:
    print("Simulation already run! Skipping...")

# Read the simulation results
print("Reading the QBlade simulation results...")

results = simulation.results(TIMESTEP_ID, BLADE_ID, C_0, RADIAL_CUTOFF, PROBE_POSITION_TOP, PROBE_POSITION_BOT, XFOIL_TRANSITION_TOP, XFOIL_TRANSITION_BOT, XFOIL_CRITICAL_AMPLIFICATION, XFOIL_MAX_ITER)

# Determine the atmospheric conditions
print("Determining the atmospheric conditions...")

turbine = simulation.turbine
z_ref = turbine.attributes["TOWERHEIGHT"]
U_ref = simulation.attributes["MEANINF"]

z_0 = surface_roughness_length(z_ref, U_ref, G, KAPPA, NU_0)
I = turbulence_intensity(z_ref, z_0)
L = turbulence_length_scale(z_ref, z_0)

# Determine the microphone positions
# TODO: AUTOMATICALLY ADD 4th row in coordinates.py function
print("Determining the microphone positions...")

turbine = simulation.turbine
z_ref = turbine.attributes["TOWERHEIGHT"]

x_t = np.array([[100.0],
                [  0.0],
                [z_ref],
                [  1.0]])

x_f = np.zeros((x_t.shape[0], x_t.shape[1], len(results["radius"])))

for i in range(x_f.shape[2]):

    matrix_fa = freestream_to_airfoil(results["alpha"][i])
    matrix_ab = airfoil_to_blade(results["radius"][i], results["chord"][i], results["twist"][i], 0.0, 0.0, 0.5, 1.0)
    matrix_bh = blade_to_hub(0.0, np.radians(0.15))
    matrix_hn = hub_to_nacelle(0.0, 0.0)
    matrix_nt = nacelle_to_turbine(turbine.attributes["TOWERHEIGHT"], 0.0, 0.0)

    matrix_ft = matrix_nt @ matrix_hn @ matrix_bh @ matrix_ab @ matrix_fa

    x_f[:, :, i] = np.linalg.inv(matrix_ft) @ x_t

x = x_f[0, :, :].T[np.newaxis, :, :]
y = x_f[1, :, :].T[np.newaxis, :, :]
z = x_f[2, :, :].T[np.newaxis, :, :]

# Determine the SPL spectra
print("Running the aeroacoustic analysis...")

f, _, _ = octave(F_MIN, F_MAX, F_REF, BASE_10)

f = f[:, np.newaxis, np.newaxis]
results["span"] = results["span"][np.newaxis, :, np.newaxis]
results["chord"] = results["chord"][np.newaxis, :, np.newaxis]
results["t_c_01"] = results["t_c_01"][np.newaxis, :, np.newaxis]
results["t_c_10"] = results["t_c_10"][np.newaxis, :, np.newaxis]
results["U"] = results["U"][np.newaxis, :, np.newaxis]
results["alpha"] = results["alpha"][np.newaxis, :, np.newaxis]
results["delta_star_top"] = results["delta_star_top"][np.newaxis, :, np.newaxis]
results["delta_star_bot"] = results["delta_star_bot"][np.newaxis, :, np.newaxis]

spl_top = te_noise(f, results["span"], results["chord"], x, y, z, results["U"], \
                   results["delta_star_top"], C_0, RHO_0, P_REF_AIR, \
                    CONVECTION_VELOCITY_COEFFICIENT, SPANWISE_CORRELATION_COEFFICIENT, \
                        True)

spl_bot = te_noise(f, results["span"], results["chord"], x, y, z, results["U"], \
                   results["delta_star_bot"], C_0, RHO_0, P_REF_AIR, \
                    CONVECTION_VELOCITY_COEFFICIENT, SPANWISE_CORRELATION_COEFFICIENT, \
                        True)

spl = 10 * np.log10(np.pow(10, spl_top/10) + np.pow(10, spl_bot/10))

spl = 10 * np.log10(np.sum(np.pow(10, spl/10), axis=1))

plot_spectrum(f[:, 0, 0], spl[:, 0], F_MIN, F_MAX, "NM80_spectrum")
