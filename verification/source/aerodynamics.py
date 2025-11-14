"""
Author:   T. Moreira da Fonte Fonseca Teles
Email:    tmoreiradafont@tudelft.nl
Date:     2025-11-14
License:  GNU GPL 3.0
"""

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd

from source.QBlade.simulation import Simulation
from source.settings import QBLADE_SIMULATION_PATH, QBLADE_RESULTS_PATH
from verification.source.settings import BEM_DTU_PATH, BEM_NREL_PATH


# Read the QBlade simulation
simulation = Simulation(QBLADE_SIMULATION_PATH)
turbine = simulation.turbine
blade = turbine.blade

# Read the simulation results
simulation.read_results(QBLADE_RESULTS_PATH)

# Read the verification data
DTU = pd.read_csv(BEM_DTU_PATH)
NREL = pd.read_csv(BEM_NREL_PATH)

# Determine the meshgrid of the radius and azimuth
azimuth = np.linspace(0.0, 2.0 * np.pi, 360)
radius = np.linspace(blade.radius[0], blade.radius[-1], 1000)
RADIUS, AZIMUTH = np.meshgrid(radius, azimuth)

# Show the steady-state operating conditions
inflow_velocity = simulation.get_results("inflow_velocity", azimuth)
tip_speed_ratio = turbine.get_results("tip_speed_ratio", azimuth)
blade_pitch = turbine.get_results("pitch", azimuth)

inflow_velocity = np.mean(inflow_velocity)
tip_speed_ratio = np.mean(tip_speed_ratio)
blade_pitch = np.mean(blade_pitch)

print(f"Inflow Velocity: {inflow_velocity} [m/s]")
print(f"Tip Speed Ratio: {tip_speed_ratio} [-]")
print(f"Blade Pitch Angle: {np.degrees(blade_pitch)} [deg]")

# Show the steady-state rotor performance
thrust_coefficient = turbine.get_results("thrust_coefficient", azimuth)
torque_coefficient = turbine.get_results("torque_coefficient", azimuth)
power_coefficient = turbine.get_results("power_coefficient", azimuth)

thrust_coefficient = np.mean(thrust_coefficient)
torque_coefficient = np.mean(torque_coefficient)
power_coefficient = np.mean(power_coefficient)

print(f"Thrust coefficient: {thrust_coefficient} [-]")
print(f"Torque coefficient: {torque_coefficient} [-]")
print(f"Power coefficient: {power_coefficient} [-]")

# Show the steady-state force distributions
axial_force = blade.get_results("axial_force", AZIMUTH, RADIUS)
tangential_force = blade.get_results("tangential_force", AZIMUTH, RADIUS)

axial_force = np.mean(axial_force, axis=0)
tangential_force = np.mean(tangential_force, axis=0)

plt.plot(radius, axial_force, label="QBlade")
plt.plot(DTU["radius"], DTU["F_x"], ls="--", label="HAWC2")
plt.plot(NREL["radius"], NREL["F_x"], ls="--", label="OpenFAST")
plt.xlabel("Radius, [m]")
plt.ylabel("Axial Force per Unit Span, [N/m]")
plt.xlim(blade.radius[0], blade.radius[-1])
plt.ylim(0.0, 12000.0)
plt.grid()
plt.legend()
plt.show()

plt.plot(radius, tangential_force, label="QBlade")
plt.plot(DTU["radius"], -DTU["F_y"], ls="--", label="HAWC2")
plt.plot(NREL["radius"], -NREL["F_y"], ls="--", label="OpenFAST")
plt.xlabel("Radius, [m]")
plt.ylabel("Tangential Force per Unit Span, [N/m]")
plt.xlim(blade.radius[0], blade.radius[-1])
plt.ylim(0.0, 1400.0)
plt.grid()
plt.legend()
plt.show()

# Show the steady-state blade deflections
axial_deflection = blade.get_results("axial_deflection", AZIMUTH, RADIUS)
radial_twist = blade.get_results("radial_twist", AZIMUTH, RADIUS)

axial_deflection = np.mean(axial_deflection, axis=0)
radial_twist = np.mean(radial_twist, axis=0)

plt.plot(radius, axial_deflection, label="QBlade")
plt.plot(DTU["radius"], DTU["delta_x"], ls="--", label="HAWC2")
plt.plot(NREL["radius"], NREL["delta_x"], ls="--", label="OpenFAST")
plt.xlabel("Radius, [m]")
plt.ylabel("Blade Deflection, [m]")
plt.xlim(blade.radius[0], blade.radius[-1])
plt.ylim(0.0, 20.0)
plt.grid()
plt.legend()
plt.show()

plt.plot(radius, np.degrees(radial_twist), label="QBlade")
plt.plot(DTU["radius"], DTU["theta_z"], ls="--", label="HAWC2")
plt.plot(NREL["radius"], NREL["theta_z"], ls="--", label="OpenFAST")
plt.xlabel("Radius, [m]")
plt.ylabel(r"Blade Twist, [$^\circ$]")
plt.xlim(blade.radius[0], blade.radius[-1])
plt.ylim(-6.0, 0.0)
plt.grid()
plt.legend()
plt.show()
