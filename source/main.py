"""
Author:   T. Moreira da Fonte Fonseca Teles
Email:    tmoreiradafont@tudelft.nl
Date:     2025-11-14
License:  GNU GPL 3.0
"""

import os

import matplotlib.pyplot as plt
import numpy as np

from source.acoustics.frequency import one_third_octave, source_frequency
from source.acoustics.leading_edge_noise import leading_edge_noise
from source.acoustics.trailing_edge_noise import trailing_edge_noise
from source.atmosphere.absorption import attenuation_coefficient, molar_concentration
from source.atmosphere.ISA import density, kinematic_viscosity, speed_of_sound
from source.atmosphere.turbulence import turbulence_intensity, turbulence_length_scale, surface_roughness_length
from source.boundary_layer.displacement_thickness import displacement_thickness
from source.coordinates.retarded_coordinates import retarded_coordinates
from source.QBlade.qblade import QBlade
from source.QBlade.simulation import Simulation
from source.settings import QBLADE_DLL_PATH, QBLADE_OPENCL_DEVICE, QBLADE_OPENCL_GROUP_SIZE
from source.settings import QBLADE_SIMULATION_PATH, QBLADE_RESULTS_PATH
from source.settings import BASE_10, F_MIN, F_MAX
from source.settings import ASPECT_RATIO, N_AZIMUTH, RADIUS_CUTOFF
from source.settings import P_0, T_0, H_R
from source.settings import XFOIL_EXE_PATH, XFOIL_ITERATION_LIMIT, XFOIL_PROBE_UPPER, XFOIL_PROBE_LOWER
from source.settings import XFOIL_TRANSITION_UPPER, XFOIL_TRANSITION_LOWER, XFOIL_AMPLIFICATION
from source.settings import N_OBSERVER, R_OBSERVER, Z_OBSERVER


# Run the QBlade simulation
print("Running the QBlade simulation...")

if not os.path.exists(QBLADE_RESULTS_PATH):

    # Create the results directory
    results_directory = os.path.dirname(QBLADE_RESULTS_PATH)
    results_name = os.path.splitext(os.path.basename(QBLADE_RESULTS_PATH))[0]
    os.makedirs(results_directory)

    # Initialise the simulation
    qblade = QBlade(QBLADE_DLL_PATH)
    qblade.create_instance(QBLADE_OPENCL_DEVICE, QBLADE_OPENCL_GROUP_SIZE)
    qblade.load_sim_definition(QBLADE_SIMULATION_PATH)
    qblade.initialise_simulation()

    # Run the simulation
    qblade.run_full_simulation()
    qblade.export_results(0, results_directory, results_name, "")

    # Close the simulation
    qblade.close_instance()

# Run the aeroacoustic simulation
print("Running the aeroacoustic simulation...")

if True:

    # Read the simulation data
    simulation = Simulation(QBLADE_SIMULATION_PATH)
    turbine = simulation.turbine
    blade = turbine.blade

    # Read the simulation results
    simulation.read_results(QBLADE_RESULTS_PATH)

    # Determine the one-third octave frequencies
    f_o, _, _ = one_third_octave(F_MIN, F_MAX, BASE_10)

    # Determine the rotor blade discretisation
    r, b, c = blade.discretise(ASPECT_RATIO, RADIUS_CUTOFF)

    # Determine the blade thicknesses
    t_01 = blade.get_geometry("thickness_01", r)
    t_10 = blade.get_geometry("thickness_10", r)

    # Determine the observer coordinates
    x_o_t = np.array([[np.cos(np.linspace(0.0, 2 * np.pi, N_OBSERVER + 1)) * R_OBSERVER],
                      [np.sin(np.linspace(0.0, 2 * np.pi, N_OBSERVER + 1)) * R_OBSERVER],
                      [np.ones(N_OBSERVER + 1) * Z_OBSERVER]])

    # Determine the atmospheric conditions
    rho_0 = density(P_0, T_0)
    nu_0 = kinematic_viscosity(T_0, rho_0)
    c_0 = speed_of_sound(T_0)

    # Determine the attenuation coefficient
    h_0 = molar_concentration(P_0, T_0, H_R)
    alpha_0 = attenuation_coefficient(f_o, P_0, T_0, h_0)

    # Determine the azimuthal and radial discretisations
    psi = np.linspace(0.0, 2 * np.pi, N_AZIMUTH)
    R, PSI = np.meshgrid(r, psi)

    # Determine the azimuthally-averaged distributions
    U = np.mean(blade.get_results("velocity", PSI, R), axis=0)
    alpha = np.mean(blade.get_results("angle_of_attack", PSI, R), axis=0)

    # Determine the displacement thicknesses
    delta_star_upper, delta_star_lower = displacement_thickness(blade, c, r, U, alpha, nu_0, XFOIL_EXE_PATH, \
                                         XFOIL_ITERATION_LIMIT, XFOIL_TRANSITION_UPPER, XFOIL_TRANSITION_LOWER, \
                                         XFOIL_AMPLIFICATION, XFOIL_PROBE_UPPER, XFOIL_PROBE_LOWER)

    # Initialise the SPL arrays
    spl_le = np.empty((N_AZIMUTH, len(f_o), len(r), N_OBSERVER + 1))
    spl_te_upper = np.empty((N_AZIMUTH, len(f_o), len(r), N_OBSERVER + 1))
    spl_te_lower = np.empty((N_AZIMUTH, len(f_o), len(r), N_OBSERVER + 1))

    # Iterate over the azimuth angles
    for i in range(N_AZIMUTH):

        # Determine the azimuth angle
        psi = 2.0 * np.pi * i / N_AZIMUTH

        # Determine the retarded coordinates
        x_s_t, x_o_f = retarded_coordinates(turbine, psi, r, x_o_t, c_0*100000)

        # Determine the source frequency
        f_s = source_frequency(simulation, psi, f_o, x_s_t, x_o_t, c_0)

        # Determine the turbine operating conditions
        u_hub = simulation.get_results("inflow_velocity", psi)
        z_hub = turbine.tower_height

        # Determine the turbulence characteristics
        z_0 = surface_roughness_length(u_hub, z_hub, nu_0)
        I = turbulence_intensity(x_s_t[2, :, :], z_0)
        L = turbulence_length_scale(x_s_t[2, :, :], z_0)

        # Determine the observer coordinates
        x = x_o_f[0, :, :]
        y = x_o_f[1, :, :]
        z = x_o_f[2, :, :]

        # Determine the noise SPLs
        spl_le[i, :, :, :] = leading_edge_noise(f_s, b, c, I, L, t_01, t_10, U, alpha, x, y, z, c_0, rho_0)
        spl_te_upper[i, :, :, :] = trailing_edge_noise(f_s, b, c, U, delta_star_upper, x, y, z, c_0, rho_0)
        spl_te_lower[i, :, :, :] = trailing_edge_noise(f_s, b, c, U, delta_star_lower, x, y, z, c_0, rho_0)

        # Subtract the atmospheric absorbtion
        distance = np.linalg.norm(x_o_t - x_s_t, axis=0)[np.newaxis, :, :]

        spl_le[i, :, :, :] -= alpha_0[:, np.newaxis, np.newaxis] * distance
        spl_te_upper[i, :, :, :] -= alpha_0[:, np.newaxis, np.newaxis] * distance
        spl_te_lower[i, :, :, :] -= alpha_0[:, np.newaxis, np.newaxis] * distance

    # Determine the combined SPLs
    spl_te = 10.0 * np.log10(np.pow(10.0, spl_te_upper / 10.0) + np.pow(10.0, spl_te_lower / 10.0))
    spl = 10.0 * np.log10(np.pow(10.0, spl_le / 10.0) + np.pow(10.0, spl_te / 10.0))

    # Determine the blade SPLs
    spl_blade = 10.0 * np.log10(np.sum(np.pow(10.0, spl / 10.0), axis=2))
    spl_le_blade = 10.0 * np.log10(np.sum(np.pow(10.0, spl_le / 10.0), axis=2))
    spl_te_blade = 10.0 * np.log10(np.sum(np.pow(10.0, spl_te / 10.0), axis=2))
    spl_te_upper_blade = 10.0 * np.log10(np.sum(np.pow(10.0, spl_te_upper / 10.0), axis=2))
    spl_te_lower_blade = 10.0 * np.log10(np.sum(np.pow(10.0, spl_te_lower / 10.0), axis=2))

    # Deterine the azimuthally averaged SPLs
    spl_turbine = 10.0 * np.log10(np.mean(np.pow(10.0, spl_blade / 10.0), axis=0))
    spl_le_turbine = 10.0 * np.log10(np.mean(np.pow(10.0, spl_le_blade / 10.0), axis=0))
    spl_te_turbine = 10.0 * np.log10(np.mean(np.pow(10.0, spl_te_blade / 10.0), axis=0))
    spl_te_upper_turbine = 10.0 * np.log10(np.mean(np.pow(10.0, spl_te_upper_blade / 10.0), axis=0))
    spl_te_lower_turbine = 10.0 * np.log10(np.mean(np.pow(10.0, spl_te_lower_blade / 10.0), axis=0))

    # Determine the OSPLs
    ospl = 10.0 * np.log10(np.sum(np.pow(10.0, spl_turbine / 10.0), axis=0))
    ospl_le = 10.0 * np.log10(np.sum(np.pow(10.0, spl_le_turbine / 10.0), axis=0))
    ospl_te = 10.0 * np.log10(np.sum(np.pow(10.0, spl_te_turbine / 10.0), axis=0))
    ospl_te_upper = 10.0 * np.log10(np.sum(np.pow(10.0, spl_te_upper_turbine / 10.0), axis=0))
    ospl_te_lower = 10.0 * np.log10(np.sum(np.pow(10.0, spl_te_lower_turbine / 10.0), axis=0))

    # Print the OSPLs
    print(f"Total SPL: {ospl} [dB]")
    print(f"LE SPL: {ospl_le} [dB]")
    print(f"TE SPL: {ospl_te} [dB]")
    print(f"TE SPL - Upper: {ospl_te_upper} [dB]")
    print(f"TE SPL - Lower: {ospl_te_lower} [dB]")

    # Show the spectra of the first obsever
    plt.semilogx(f_o, spl_turbine[:, 0], label="Total Noise")
    plt.semilogx(f_o, spl_le_turbine[:, 0], label="LE Noise")
    plt.semilogx(f_o, spl_te_turbine[:, 0], label="TE Noise")
    plt.semilogx(f_o, spl_te_upper_turbine[:, 0], label="TE Noise - Upper")
    plt.semilogx(f_o, spl_te_lower_turbine[:, 0], label="TE Noise - Lower")
    plt.xlabel("Frequency, [Hz]")
    plt.ylabel("SPL, [dB]")
    plt.xlim(F_MIN, F_MAX)
    plt.ylim(0.0, 60.0)
    plt.grid(which="both")
    plt.legend()
    plt.show()

    # Show the directivity pattern
    phi_observer = np.linspace(0.0, 2.0 * np.pi, x_o_t.shape[2])

    plt.polar(phi_observer, ospl, label="Total Noise")
    plt.polar(phi_observer, ospl_le, label="LE Noise")
    plt.polar(phi_observer, ospl_te, label="TE Noise")
    plt.polar(phi_observer, ospl_te_upper, label="TE Noise - Upper")
    plt.polar(phi_observer, ospl_te_lower, label="TE Noise - Lower")
    plt.ylim(0.0, 65.0)
    plt.legend()
    plt.show()
