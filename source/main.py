"""
Author:   T. Moreira da Fonte Fonseca Teles
Email:    tmoreiradafont@tudelft.nl
Date:     2025-03-18
License:  GNU GPL 3.0

Main script.

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

import numpy as np
import pandas as pd

from coordinates import turbine_to_nacelle, nacelle_to_hub, hub_to_blade, blade_to_airfoil
from inflow_noise import inflow_noise
from misc import octave
from plot import spectrum, directivity, map
from QBlade.dll import QBlade
from QBlade.simulation import Simulation
from settings import *
from trailing_edge_noise import trailing_edge_noise


# Meeting 19/02/2025
DEMO_PLOT = 2
DEMO_SOURCE = 1

# Copy files to the working directory
turbine_dir = os.path.dirname(TURBINE_PATH)
turbine_name = os.path.splitext(os.path.basename(TURBINE_PATH))[0]

simulation_name = os.path.splitext(os.path.basename(SIMULATION_PATH))[0]

turbine_work_dir = os.path.join(WORKING_DIR, simulation_name, turbine_name)
simulation_work_path = os.path.join(WORKING_DIR, simulation_name, simulation_name + ".sim")

if not os.path.exists(os.path.join(WORKING_DIR, simulation_name)):
    shutil.copytree(turbine_dir, turbine_work_dir)
    shutil.copyfile(SIMULATION_PATH, simulation_work_path)

# Configure the simulation
simulation = Simulation(simulation_work_path)
turbine = simulation.turbine

turbine.write("BEMSPEEDUP", 100.0)   # Speedup convergence
turbine.write("STRUCTURALFILE", "")  # Disable structural simulation
turbine.write("CONTROLLERFILE", "")  # Disable controller simulation
turbine.write("PARAMETERFILE", "")   # Disable controller simulation
turbine.write("CONTROLLERTYPE", 0)   # Disable controller simulation

# Run QBlade simulation
results_dir = os.path.join(WORKING_DIR, simulation_name)
results_name = simulation_name

if not os.path.exists(os.path.join(results_dir, results_name + ".txt")):

    # Start QBlade
    if os.path.exists(QBLADE_PATH):
        print(f"Using shared library file: {QBLADE_PATH}!")
    else:
        print(f"No shared library file was found in {QBLADE_PATH}!")
        sys.exit(1)

    qblade = QBlade(QBLADE_PATH)
    qblade.createInstance(QBLADE_CL_DEVICE, QBLADE_GROUP_SIZE)

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
    qblade.exportResults(0, results_dir.encode("utf-8"), results_name.encode("utf-8"), b"")

    # Unload the QBlade library
    qblade.unload_library()

results = pd.read_csv(os.path.join(results_dir, results_name + ".txt"), skiprows=2, delimiter="\t").iloc[TIMESTEP]

# Determine 1/3 band octave frequencies
f, _, _ = octave(F_MIN, F_MAX, F_REF, BASE_10)

# Determine the blade distributions
n_panels = turbine.attributes["NUMPANELS"]

U     = np.zeros(n_panels)
Re    = np.zeros(n_panels)
alpha = np.zeros(n_panels)
r     = np.zeros(n_panels)
tc_01 = np.zeros(n_panels)
tc_10 = np.zeros(n_panels)

for i in range(n_panels):
    U[i] = results[f"Total_Velocity_Blade_{BLADE}_PAN_{i}_[m/s]"]
    Re[i] = results[f"Reynolds_Number_Blade_{BLADE}_PAN_{i}_[-]"]
    alpha[i] = np.radians(results[f"Angle_of_Attack_at_0.25c_Blade_{BLADE}_PAN_{i}_[deg]"])
    r[i] = results[f"Radius_Blade_{BLADE}_PAN_{i}_[m]"]
    tc_01[i] = turbine.blade.interpolate(r[i]).thickness(0.01)
    tc_10[i] = turbine.blade.interpolate(r[i]).thickness(0.10)

radiuses = np.array(turbine.blade.data["pos"])
chords = np.array(turbine.blade.data["chord"])
twists    = np.array(turbine.blade.data["twist"])
offsets_x = np.array(turbine.blade.data["offset_x"])
offsets_y = np.array(turbine.blade.data["offset_y"])
p_axes   = np.array(turbine.blade.data["p_axis"])

if turbine.attributes["DISCTYPE"] == 0:
    spans = np.diff(radiuses)
    span = np.interp(r, radiuses, spans)

elif turbine.attributes["DISCTYPE"] == 1:
    span = np.ones(n_panels) * (radiuses[-1] - radiuses[0]) / n_panels

elif turbine.attributes["DISCTYPE"] == 2:
    r_virtual = np.concat((np.array([radiuses[0]]), r, np.array([radiuses[-1]])))
    span = (r_virtual[2:n_panels+2] - r_virtual[0:n_panels]) / 2


pos = r
chord = np.interp(r, radiuses, chords)
twist = np.interp(r, radiuses, twists)
offset_x = np.interp(r, radiuses, offsets_x)
offset_y = np.interp(r, radiuses, offsets_y)
p_axis = np.interp(r, radiuses, p_axes)

if DEMO_PLOT == 0:
    x_t = np.array([[100],
                    [0],
                    [turbine.attributes["TOWERHEIGHT"]]])

elif DEMO_PLOT == 1:

    N = 73
    NN = 36

    theta = np.linspace(0, 2*np.pi, N)[:, np.newaxis]
    azi = np.linspace(0, 2*np.pi, NN, endpoint=False)[np.newaxis, :]

    x = 100 * np.cos(theta)
    y = 100 * np.sin(theta)
    z = 1.2 * np.ones(N)[:, np.newaxis]
    
    x = x
    y = y
    z = z - 178.4

    x_1 = (x * np.ones(azi.shape))
    y_1 = np.cos(azi) * y - np.sin(azi) * z
    z_1 = np.sin(azi) * y + np.cos(azi) * z

    x_2 = x_1
    y_2 = y_1
    z_2 = z_1 + 178.4

    x_t = np.array([x_2.flatten(),
                    y_2.flatten(),
                    z_2.flatten()])

elif DEMO_PLOT == 2:

    N = 37

    theta = np.linspace(0, 2*np.pi, N)

    locs_y = turbine.attributes["TOWERHEIGHT"] * np.cos(theta)
    locs_z = turbine.attributes["TOWERHEIGHT"] * (1 - np.sin(theta))

    x_t = np.array([20 * np.ones(N),
                    locs_y,
                    locs_z])

else:
    print("Unknown demonstration!")
    sys.exit(1)

x_t = x_t[np.newaxis, :]

x_n = turbine_to_nacelle(x_t, turbine.attributes["TOWERHEIGHT"], turbine.attributes["SHAFTTILT"], 0)
x_h = nacelle_to_hub(x_n, turbine.attributes["OVERHANG"], 0)
x_b = hub_to_blade(x_h, turbine.attributes["ROTORCONE"], 0)
x_a = blade_to_airfoil(x_b, pos, chord, twist, offset_x, offset_y, p_axis)

x = x_a[:, 0, :]
y = x_a[:, 1, :]
z = x_a[:, 2, :]

# Add np.newaxis to match the dimensions
f = f[:, np.newaxis, np.newaxis]
span = span[np.newaxis, :, np.newaxis]
chord = chord[np.newaxis, :, np.newaxis]
tc_01 = tc_01[np.newaxis, :, np.newaxis]
tc_10 = tc_10[np.newaxis, :, np.newaxis]
pos = pos[np.newaxis, :, np.newaxis]
x = x[np.newaxis, :, :]
y = y[np.newaxis, :, :]
z = z[np.newaxis, :, :]
U = U[np.newaxis, :, np.newaxis]
Re = Re[np.newaxis, :, np.newaxis]
alpha = alpha[np.newaxis, :, np.newaxis]

# Determine the SPL spectra
if DEMO_SOURCE == 0:
    z_g = 178.4 * np.ones(n_panels)[np.newaxis, :, np.newaxis] # TODO: Use the correct value
    spl = inflow_noise(f, span, chord, tc_01, tc_10, x, y, z, U, alpha, C_0, RHO_0, z_g, Z_0, True, "ZHS", "ZHS")

elif DEMO_SOURCE == 1:
    spl = trailing_edge_noise(turbine.blade, f, pos, x, y, z, U, Re, alpha, span, chord, C_0, RHO_0, P_REF, PROBE_TOP, PROBE_BOT, RADIAL_CUTOFF, MAX_ITER)

if DEMO_PLOT == 0:
    spectrum(f.squeeze(), spl, 3, 100)

elif DEMO_PLOT == 1:
    directivity(spl, 3, NN, N, 1.2, 100)

elif DEMO_PLOT == 2:
    map(spl, n_panels, N, 20)
