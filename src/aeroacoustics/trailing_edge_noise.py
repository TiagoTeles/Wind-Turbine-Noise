"""
Author:   T. Moreira da Fonte Fonseca Teles
Email:    tmoreiradafont@tudelft.nl
Date:     2025-02-19
License:  GNU GPL 3.0

Calculate the TBLTE noise.

Classes:
    None

Functions:
    TBD

Exceptions:
    None
"""

import os
import sys

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import scipy as sp

sys.path.append(os.path.dirname(sys.path[0]))

from octave import octave
from QBlade.turbine import Turbine
from XFoil.xfoil import XFoil


BASE_10 = True          # Use base 10?
F_MIN = 20              # Minimum frequency, [Hz]
F_MAX = 20000           # Maximum frequency, [Hz]
F_REF = 1000            # Reference frequency, [Hz]
P_REF = 2E-5            # Reference pressure, [Pa]
C_0 = 340.0             # Speed of sound, [m/s]        TODO: Base on temperature?
RHO_0 = 1.225           # Air density, [kg/m^3]        TODO: Use QBalde value?
BLADE = 1               # Blade id                     TODO: Does not work for unsteady case
TIMESTEP = -1           # Timestep index to use
XFOIL_PATH = "bin\\XFoil\\xfoil.exe"
PROBE_TOP = 0.975
PROBE_BOT = 0.950


def E(x):

    S_2, C_2 = sp.special.fresnel(np.sqrt(2 * x / np.pi))

    return C_2 - 1j * S_2


def roger_moreau(blade, f, x, y, z, r, U, Re, alpha, b, c, c_0, rho_0, P_REF, xfoil_path):

    
    # Determine the number of frequencies and panels
    n_panels = U.size


    delta_f = 0.232 * f

    omega = 2 * np.pi * f

    M = U / C_0

    beta = np.sqrt(1 - np.square(M))

    S_0 = np.sqrt(np.square(x) + np.square(beta) * (np.square(y) + np.square(z)))

    K_x = omega / U

    delta_star_top = np.zeros(n_panels)
    delta_star_bot = np.zeros(n_panels)

    # Iterate through each panel
    for i in range(n_panels):

        # Determine the neighbouring airfoils
        index_inboard  = blade.data[blade.data["pos"] < r[i]].index.max()
        index_outboard = blade.data[blade.data["pos"] > r[i]].index.min()

        r_inboard  = blade.data["pos"].iloc[index_inboard]
        r_outboard = blade.data["pos"].iloc[index_outboard]

        polar_path_inboard  = blade.data["polar_path"].iloc[index_inboard]
        polar_path_outboard = blade.data["polar_path"].iloc[index_outboard]

        airfoil_path_inboard  = None
        airfoil_path_outboard = None

        for polar in blade.polars:

            if polar.path == polar_path_inboard:
                airfoil_path_inboard = polar.attributes["FOILNAME"]

            if polar.path == polar_path_outboard:
                airfoil_path_outboard = polar.attributes["FOILNAME"]

        if airfoil_path_inboard is None or airfoil_path_outboard is None:
            print("Polar not found!")
            sys.exit(1)

        # TODO: PRECOMPUTE THE INTERPOLATED AIRFOILS
        # TODO: INTERPOLATE BETWEEN THE INBOARD AND OUTBOARD AIRFOILS TO GET THE CORRECT AIRFOIL
        # TODO: Is suction side always the top side? Is pressure side always the bottom side?
        if "Circular" not in airfoil_path_inboard:
            xfoil = XFoil(xfoil_path)
            top, bot = xfoil.run(airfoil_path_inboard, Re[i], M[i], alpha[i], it=100)

            delta_star_top[i] = np.interp(PROBE_TOP, top["x/c"], top["delta_star"])
            delta_star_bot[i] = np.interp(PROBE_BOT, bot["x/c"], bot["delta_star"])

        else:
            delta_star_top[i] = 1E-9
            delta_star_bot[i] = 1E-9

    # TODO: For now average the delta_star_top and delta_star_bot values, and use factor of 4 in SPL calculation
    # TODO: This is because for a flat plate delta_star_top == delta_star_bot

    delta_star = (delta_star_top + delta_star_bot) / 2

    omega_tilde = K_x * delta_star

    # TODO: IS THIS THE BEST WAY TO CALCULATE F? IS THERE A NEWER MODEL
    F = (33.28 * omega_tilde) / (1 - 5.489 * omega_tilde + 36.74 * np.square(omega_tilde) + 0.1505 * np.pow(omega_tilde, 5))

    # TODO: ANY NEWER EMPITICAL FORMULATION? IS 2E-5 related to P_REF???
    Phi_pp = np.square(0.5 * rho_0 * np.square(U)) * (delta_star / U) * 2E-5 * F

    # TODO: FIND VALUE OF b_c. Roger claims it is a constant, but does not provide a value
    b_c = 1

    # TODO: FIND VALUE OF alpha_c. Roger does not provide a value
    alpha_c = 1 / 0.8

    U_c = U / alpha_c

    k = omega / c_0

    K_2 = k * y / S_0

    l_y = (omega / (b_c * U_c)) / (np.square(K_2) + np.square(omega) / np.square(b_c * U_c))

    # Prof. Casalino correction term
    l_y *= (1 - np.exp(-np.square(omega_tilde) / 0.09))

    Pi_0 = (1/np.pi) * Phi_pp * l_y

    # TODO: Why not U_c
    K = omega / U

    K_line = K * (c/2)

    mu_line = K_line * M / np.square(beta)

    epsilon = np.pow(1 + 1 / (4 * mu_line), -0.5)

    K_2_line = K_2 * (c/2)

    kappa_line = np.sqrt(np.square(mu_line) - np.square(K_2_line) / np.square(beta))

    if np.any(np.square(mu_line) - np.square(K_2_line) / np.square(beta) < 0):
        print("Sub-critical gust detected!")
        print("I probably need to implement them :(")


    K_1_line = K_line * alpha_c

    B = K_1_line + M * mu_line + kappa_line
    C = K_1_line - mu_line * (x / S_0 - M)

    # Determine supercritical gust contribution
    # Remove  "-np.exp(-2*1j*C)" as indicated in paper
    # TODO: DOES THIS PART NEED TO BE CORRECTED? HAS {} but not ^c
    uncorrected_1 = (1 + 1j) * np.exp(-2 * 1j * C) * np.sqrt(B/(B-C)) * E(2 * (B - C)) - (1 + 1j) * E(2*B) +1

    corrected_1 = np.zeros(uncorrected_1.shape, dtype=complex)
    corrected_1.real = uncorrected_1.real
    corrected_1.imag = uncorrected_1.imag * epsilon

    part_1 = - np.exp(2 * 1j * C) / (1j * C) * corrected_1

    uncorrected_2 = np.exp(4 * 1j * kappa_line) * (1 - (1 + 1j) * E(4 * kappa_line))

    corrected_2 = np.zeros(uncorrected_2.shape, dtype=complex)
    corrected_2.real = uncorrected_2.real
    corrected_2.imag = uncorrected_2.imag * epsilon

    D = kappa_line - mu_line * x / S_0

    G = (1 + epsilon) * np.exp(1j * (2 * kappa_line + D)) * np.sin(D - 2 * kappa_line) / (D - 2 * kappa_line) + (1 - epsilon) * np.exp(1j * (-2 * kappa_line + D)) * np.sin(D + 2*kappa_line) / (D + 2 * kappa_line)
    G += (1 + epsilon) * (1 - 1j) / (2 * (D - 2 * kappa_line)) * np.exp(4 * 1j * kappa_line) * E(4 * kappa_line) - (1 - epsilon) * (1 + 1j) / (2 * (D + 2 * kappa_line)) * np.exp(-4 * 1j * kappa_line) * E(4 * kappa_line)
    G += np.exp(2 * 1j * D) / 2 * np.sqrt(2 * kappa_line / D) * E(2 * D) * ((1 + 1j) * (1 - epsilon) / (D + 2 * kappa_line) - (1 - 1j) * (1 + epsilon) / (D - 2 * kappa_line))

    part_2 = corrected_2 - np.exp(2 * 1j * D) + 1j * (D + K_line + M * mu_line - kappa_line) * G

    Theta = np.sqrt((K_1_line + M * mu_line + kappa_line) / (K_line + M * mu_line + kappa_line))

    H = (1 + 1j) * np.exp(-4 * 1j * kappa_line) * (1 - np.square(Theta)) / (2 * np.sqrt(np.pi) * (alpha_c - 1) * K_line * np.sqrt(B))

    I = part_1 + H * part_2

    S_pp = np.square(omega * z * c / (4 * np.pi * c_0 * np.square(S_0))) \
           * 2 * np.pi * b * np.square(np.abs(I)) * Pi_0

    SPL = 10 * np.log10(4 * np.pi * delta_f * S_pp / np.square(P_REF))

    return SPL

# TODO: COMPUTE SUBCRITICAL GUST AND APPLY CONTINUITY CONDITION AT CRITICAL GUST
# TODO: IMPLEMENT CORRECTIONS FROM PAPER FROM PROF DAMIANO CASALINO
# TODO: DETERMINE CORRECT FRESENL INTEGRAL
if __name__ == "__main__":

    # Load data
    turbine = Turbine("tmp\\QBlade\\DTU_10MW\\DTU_10MW_RWT.trb")
    results = pd.read_csv("tmp\\QBlade\\DTU_10MW_RWT.txt", skiprows=2, delimiter="\t").iloc[TIMESTEP]

    # Determine 1/3 band octave frequencies
    f, _, _ = octave(F_MIN, F_MAX, F_REF, BASE_10)

    f = f[:, np.newaxis]

    # Determine the number of panels
    n_panels = turbine.attributes["NUMPANELS"]

    # Determine the blade properties
    chords = np.array(turbine.blade.data["chord"])
    radiuses = np.array(turbine.blade.data["pos"])

    # Initialise inputs
    U     = np.zeros(n_panels)
    Re    = np.zeros(n_panels)
    alpha = np.zeros(n_panels)
    r     = np.zeros(n_panels)

    for i in range(n_panels):
        U[i] = results[f"Total_Velocity_Blade_1_PAN_{i}_[m/s]"]
        Re[i] = results[f"Reynolds_Number_Blade_{BLADE}_PAN_{i}_[-]"]
        alpha[i] = np.radians(results[f"Angle_of_Attack_at_0.25c_Blade_{BLADE}_PAN_{i}_[deg]"])
        r[i] = results[f"Radius_Blade_{BLADE}_PAN_{i}_[m]"]

    if turbine.attributes["DISCTYPE"] == 0:
        b = np.diff(radiuses)

    elif turbine.attributes["DISCTYPE"] == 1:
        b = np.ones(n_panels) * (radiuses[-1] - radiuses[0]) / n_panels

    elif turbine.attributes["DISCTYPE"] == 2:
        r_virtual = np.concat((np.array([radiuses[0]]), r, np.array([radiuses[-1]])))
        b = (r_virtual[2:n_panels+2] - r_virtual[0:n_panels]) / 2

    else:
        print("Unkown discretisation type!")
        sys.exit(1)

    c = np.interp(r, radiuses, chords)

    spl = roger_moreau(turbine.blade, f, 0.0, 0.0, 1.0, r, U, Re, alpha, b, c, C_0, RHO_0, P_REF, XFOIL_PATH)

    # Plot spectra
    for i in range(spl.shape[1]):
        plt.plot(f, spl[:, i], label=f"Panel {i}")

    total_ten = 10 * np.log10(np.sum(np.pow(10, spl/10), axis=1))

    plt.plot(f, total_ten, label="Total TEN", lw=2)

    total_inflow = np.array([86.26640703,   84.94987456,
                             83.52741888,   81.99146557,   80.33292237,   78.54085755,   76.60207578,
                             74.50064419,   72.21743921,   69.72979754,   67.01135896,   64.03218189,
                             60.75917678,   57.15678214,   53.18748525,   48.81113156,   43.98110353,
                             38.63524134,   32.68118276,   25.97931803,   18.32797496,    9.45156256,
                             -1.01284176,  -13.53454012,  -28.71239827,  -47.29850333,  -70.22503201,
                            -98.64437206, -133.99007291, -178.05819394, -233.10720716])

    plt.plot(f, total_inflow, label="Total INFLOW", lw=2)

    plt.plot(f, 10 * np.log10(np.pow(10, total_ten / 10) + np.pow(10, total_inflow / 10)), label="Total", lw=2)

    plt.xlim(F_MIN, F_MAX)
    plt.ylim(0, 100)
    plt.xlabel("Frequency, [Hz]")
    plt.ylabel("SPL, [dB]")
    plt.xscale("log")
    plt.yscale("linear")
    plt.legend()
    plt.grid(which="both")
    plt.show()
