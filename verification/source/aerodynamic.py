"""
Author:   T. Moreira da Fonte Fonseca Teles
Email:    tmoreiradafont@tudelft.nl
Date:     2025-10-20
License:  GNU GPL 3.0
"""

import pandas as pd
import matplotlib.pyplot as plt
import numpy as np

from source.QBlade.simulation import Simulation
from source.settings import QBLADE_SIMULATION_PATH, QBLADE_RESULTS_PATH


# Read the simulation results
simulation = Simulation(QBLADE_SIMULATION_PATH)
simulation.read_results(QBLADE_RESULTS_PATH)

# Read the verification data
data_hawc2 = pd.read_csv("verification/data/aerodynamics/HAWC2.csv")
data_openfast = pd.read_csv("verification/data/aerodynamics/OpenFAST.csv")

# Determine the azimuth angles
azimuth = np.linspace(0.0, 2.0 * np.pi, 180)

# Determine the blade radius
root_radius = simulation.turbine.blade.radius[0]
tip_radius = simulation.turbine.blade.radius[-1]
radius = np.linspace(root_radius, tip_radius, 1000)

# Show the steady-state operating conditions
inflow_velocity = simulation.get_results("inflow_velocity", azimuth)
tip_speed_ratio = simulation.turbine.get_results("tip_speed_ratio", azimuth)
blade_pitch = simulation.turbine.get_results("pitch", azimuth)

print(f"Inflow Velocity: {np.mean(inflow_velocity)} [m/s]")
print(f"Tip Speed Ratio: {np.mean(tip_speed_ratio)} [-]")
print(f"Blade Pitch Angle: {np.degrees(np.mean(blade_pitch))} [deg]")

# Show the steady-state rotor performance
thrust_coefficient = simulation.turbine.get_results("thrust_coefficient", azimuth)
torque_coefficient = simulation.turbine.get_results("torque_coefficient", azimuth)
power_coefficient = simulation.turbine.get_results("power_coefficient", azimuth)

print(f"Thrust coefficient: {np.mean(thrust_coefficient)} [-]")
print(f"Torque coefficient: {np.mean(torque_coefficient)} [-]")
print(f"Power coefficient: {np.mean(power_coefficient)} [-]")

# Show the steady-state force distribution
axial_force = simulation.turbine.blade.get_results("axial_force", azimuth, radius)
tangential_force = simulation.turbine.blade.get_results("tangential_force", azimuth, radius)

plt.plot(radius, np.mean(axial_force, axis=0), label="QBlade")
plt.plot(data_hawc2["radius"], data_hawc2["F_x"], ls="--", label="HAWC2")
plt.plot(data_openfast["radius"], data_openfast["F_x"], ls="--", label="OpenFAST")
plt.xlabel("Radius, [m]")
plt.ylabel("Axial Force per Unit Span, [N/m]")
plt.xlim(root_radius, tip_radius)
plt.ylim(0.0, 12000.0)
plt.grid()
plt.legend()
plt.show()

plt.plot(radius, np.mean(tangential_force, axis=0), label="QBlade")
plt.plot(data_hawc2["radius"], -data_hawc2["F_y"], ls="--", label="HAWC2")
plt.plot(data_openfast["radius"], -data_openfast["F_y"], ls="--", label="OpenFAST")
plt.xlabel("Radius, [m]")
plt.ylabel("Tangential Force per Unit Span, [N/m]")
plt.xlim(root_radius, tip_radius)
plt.ylim(0.0, 1400.0)
plt.grid()
plt.legend()
plt.show()

# Show the steady-state blade deflections
axial_deflection = simulation.turbine.blade.get_results("axial_deflection", azimuth, radius)
radial_twist = simulation.turbine.blade.get_results("radial_twist", azimuth, radius)

plt.plot(radius, np.mean(axial_deflection, axis=0), label="QBlade")
plt.plot(data_hawc2["radius"], data_hawc2["delta_x"], ls="--", label="HAWC2")
plt.plot(data_openfast["radius"], data_openfast["delta_x"], ls="--", label="OpenFAST")
plt.xlabel("Radius, [m]")
plt.ylabel("Blade Deflection, [m]")
plt.xlim(root_radius, tip_radius)
plt.ylim(0.0, 20.0)
plt.grid()
plt.legend()
plt.show()

plt.plot(radius, np.degrees(np.mean(radial_twist, axis=0)), label="QBlade")
plt.plot(data_hawc2["radius"], data_hawc2["theta_z"], ls="--", label="HAWC2")
plt.plot(data_openfast["radius"], data_openfast["theta_z"], ls="--", label="OpenFAST")
plt.xlabel("Radius, [m]")
plt.ylabel(r"Blade Twist, [$^\circ$]")
plt.xlim(root_radius, tip_radius)
plt.ylim(-6.0, 0.0)
plt.grid()
plt.legend()
plt.show()
