"""
Author:   T. Moreira da Fonte Fonseca Teles
Email:    tmoreiradafont@tudelft.nl
Date:     2025-03-21
License:  GNU GPL 3.0
"""

import os
import shutil

import numpy as np

from aeroacoustics.inflow_noise import inflow_noise
from coordinates import turbine_to_nacelle, nacelle_to_hub, hub_to_blade, blade_to_airfoil
from misc import octave, surface_roughness, turbulence_intensity, turbulence_length
from plot import spectrum, directivity, map
from QBlade.simulation import Simulation
from settings import *


# Temporary constants for testing
X_SPECTRA = 100
R_DIRECTIVITY = 100
Z_DIRECTIVITY = 1.2
X_BLADE = 100
Z_BLADE = 0.0
N_MICROPHONES = 13
N_AZIMUTH = 180

# Create a temporary copy of the simulation
print("Creating a temporary copy of the simulation...")

source_dir = os.path.dirname(SIMULATION_PATH)
destination_dir = os.path.join("temp", os.path.basename(source_dir))

if not os.path.exists(destination_dir):
    shutil.copytree(source_dir, destination_dir)

# Create a simulation object
print("Reading the simulation and turbine definitions...")

simulation_name = os.path.basename(source_dir) + ".sim"
simulation_path = os.path.join(destination_dir, simulation_name)

simulation = Simulation(simulation_path)

# Configure the simulation
print("Configuring the simulation...")

turbine = simulation.turbine

turbine.write("BEMSPEEDUP", 60.0)    # Set BEM convergence acceleration time
turbine.write("STRUCTURALFILE", "")  # Disable structural simulation
turbine.write("CONTROLLERFILE", "")  # Disable controller simulation
turbine.write("PARAMETERFILE", "")   # Disable controller simulation
turbine.write("CONTROLLERTYPE", 0)   # Disable controller simulation

# Run the simulation
print("Running the simulation...")

results_name = os.path.basename(source_dir) + ".txt"
results_path = os.path.join(destination_dir, results_name)

if not os.path.exists(results_path):

    # Initialise the simulation
    simulation.initialise(QBLADE_PATH, QBLADE_CL_DEVICE, QBLADE_GROUP_SIZE)

    # Run the simulation
    simulation.run(N_TIMESTEP)

    # Close the simulation
    simulation.close()

# Read the simulation results
results = simulation.results(TIMESTEP, BLADE)

# Determine the SPL spectra
print("Running the aeroacoustic analysis...")
f, _, _ = octave(F_MIN, F_MAX, F_REF, BASE_10)

z_0 = surface_roughness(turbine.attributes["TOWERHEIGHT"], simulation.attributes["MEANINF"], G, \
                        KAPPA, NU_0)
I = turbulence_intensity(turbine.attributes["TOWERHEIGHT"], z_0)     # TODO: DETERMINE FOR EACH PANEL?
L = turbulence_length(turbine.attributes["TOWERHEIGHT"], z_0)

f = f[:, np.newaxis, np.newaxis]
results["span"] = results["span"][np.newaxis, :, np.newaxis]
results["chord"] = results["chord"][np.newaxis, :, np.newaxis]
results["tc_01"] = results["tc_01"][np.newaxis, :, np.newaxis]
results["tc_10"] = results["tc_10"][np.newaxis, :, np.newaxis]
results["U"] = results["U"][np.newaxis, :, np.newaxis]
results["aoa"] = results["aoa"][np.newaxis, :, np.newaxis]

x_n_spectra = np.array([[X_SPECTRA],
                        [0],
                        [0]])

x_h_spectra = nacelle_to_hub(x_n_spectra, turbine.attributes["OVERHANG"], 0)
x_b_spectra = hub_to_blade(x_h_spectra, turbine.attributes["ROTORCONE"], results["pitch"])
x_a_spectra = blade_to_airfoil(x_b_spectra, results["pos"], results["chord"], results["twist"], \
                               results["offset_x"], results["offset_y"], results["p_axis"])

x = x_a_spectra[:, 0, :][np.newaxis, :, :]
y = x_a_spectra[:, 1, :][np.newaxis, :, :]
z = x_a_spectra[:, 2, :][np.newaxis, :, :]

spl_inflow_spectra = inflow_noise(f, results["span"], results["chord"], results["tc_01"], \
                                  results["tc_10"], x, y, z, results["U"], results["aoa"], \
                                  I, L, C_0, RHO_0, True)

spectrum(f.squeeze(), spl_inflow_spectra, F_MIN, F_MAX, turbine.attributes["NUMBLADES"], X_SPECTRA)








yaw = np.linspace(0, 2*np.pi, N_MICROPHONES)
azimuth = np.linspace(0, 2*np.pi, N_AZIMUTH, endpoint=False)

x_t_directivity = np.array([R_DIRECTIVITY * np.cos(yaw),
                            R_DIRECTIVITY * np.sin(yaw),
                            0 * np.ones(N_MICROPHONES)])

x_n_directivity = turbine_to_nacelle(x_t_directivity, turbine.attributes["TOWERHEIGHT"], \
                               turbine.attributes["SHAFTTILT"], results["yaw"])

# CAUSE OF ASSYMETRY IN DIRECTIVITY PATTERN
# x_h_directivity = nacelle_to_hub(x_n_directivity, np.ones(N_AZIMUTH) * turbine.attributes["OVERHANG"], azimuth)
x_h_directivity = nacelle_to_hub(x_n_directivity, np.zeros(N_AZIMUTH) * turbine.attributes["OVERHANG"], azimuth)

x_h_directivity = np.transpose(np.reshape(np.transpose(x_h_directivity, (1, 0, 2)), (3, 1, N_AZIMUTH * N_MICROPHONES), order="F"), (1, 0, 2))

x_b_directivity = hub_to_blade(x_h_directivity, turbine.attributes["ROTORCONE"], results["pitch"])

x_a_directivity = blade_to_airfoil(x_b_directivity, results["pos"], results["chord"], results["twist"], \
                            results["offset_x"], results["offset_y"], results["p_axis"])

x = x_a_directivity[:, 0, :][np.newaxis, :, :]
y = x_a_directivity[:, 1, :][np.newaxis, :, :]
z = x_a_directivity[:, 2, :][np.newaxis, :, :]

spl_inflow_directivity = inflow_noise(f, results["span"], results["chord"], results["tc_01"], \
                                      results["tc_10"], x, y, z, results["U"], results["aoa"], \
                                      I, L, C_0, RHO_0, False)

directivity(spl_inflow_directivity, turbine.attributes["NUMBLADES"], \
             N_AZIMUTH, N_MICROPHONES, Z_DIRECTIVITY, R_DIRECTIVITY)








x_t_blade = np.array([[X_BLADE],
                      [0],
                      [Z_BLADE]])

x_n_blade = turbine_to_nacelle(x_t_blade, turbine.attributes["TOWERHEIGHT"], \
                               turbine.attributes["SHAFTTILT"], results["yaw"])


x_h_blade = nacelle_to_hub(x_n_blade, np.ones(N_AZIMUTH) * turbine.attributes["OVERHANG"], -0.5*np.pi + np.linspace(0, 2*np.pi, N_AZIMUTH, endpoint=False))

x_h_blade = np.transpose(x_h_blade, (2, 1, 0))

x_b_blade = hub_to_blade(x_h_blade, turbine.attributes["ROTORCONE"], results["pitch"])
x_a_blade = blade_to_airfoil(x_b_blade, results["pos"], results["chord"], results["twist"], \
                            results["offset_x"], results["offset_y"], results["p_axis"])

x = x_a_blade[:, 0, :][np.newaxis, :, :]
y = x_a_blade[:, 1, :][np.newaxis, :, :]
z = x_a_blade[:, 2, :][np.newaxis, :, :]

spl_inflow_blade = inflow_noise(f, results["span"], results["chord"], results["tc_01"], \
                                results["tc_10"], x, y, z, results["U"], results["aoa"], \
                                I, L, C_0, RHO_0, False)

map(spl_inflow_blade, spl_inflow_blade.shape[1], spl_inflow_blade.shape[2], X_BLADE, Z_BLADE)
