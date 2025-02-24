"""
Author:   T. Moreira da Fonte Fonseca Teles
Email:    tmoreiradafont@tudelft.nl
Date:     2025-02-21
License:  GNU GPL 3.0

Calculate the inflow noise.

Classes:
    None

Functions:
    amiet
    moriarty
    directivity
    inflow_noise

Exceptions:
    None
"""

import os
import sys

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd

sys.path.append(os.path.dirname(sys.path[0]))

from QBlade.turbine import Turbine
from aeroacoustics.octave import octave

# Settings
BASE_10 = False         # Use base 10?
F_MIN = 20              # Minimum frequency, [Hz]
F_MAX = 20000           # Maximum frequency, [Hz]
F_REF = 1000            # Reference frequency, [Hz]
SPL_CORRECTION = True   # Apply 10dB SPL correction?.
C_0 = 340.0             # Speed of sound, [m/s]        TODO: Base on temperature?
RHO_0 = 1.225           # Air density, [kg/m^3]        TODO: Use QBalde value?
I = 0.1                 # Turbulence intensity, [-]    TODO: Find model
L = 42                  # Turbulence length scale, [m] TODO: Find model
BLADE = 1               # Blade id                     TODO: Does not work for usnteady case
TIMESTEP = -1           # Timestep index to use

def amiet(f, b, c, U, alpha, r_e, L, I, c_0, rho_0):
    """
    Determine the SPL based on a flat plate.

    Arguments:
        f : np.array -- frequency array, [Hz]
        b : np.array -- span array, [m]
        c : np.array -- chord array, [m]
        U : np.array -- velocity array, [m/s]
        alpha : np.array -- angle of attack array, [rad]
        r_e : np.array -- distance array, [m]
        L : float -- turbulence length scale, [m]
        I : float -- turbulence intensity, [-]
        c_0 : float -- speed of sound, [m/s]
        rho_0 : float -- air density, [kg/m^3]

    Returns:
        spl_amiet : np.array -- SPL array, [dB]
    """

    # Determine the number of frequencies and panels
    n_frequencies = f.size
    n_panels = U.size

    # Initialise the SPL matrix
    spl_amiet = np.zeros((n_frequencies, n_panels))

   # Determine the high-frequency SPL
   # Use 78.4 instead of 58.4 due to units of rho_0 and c_0
    M = U / c_0
    K_x = (2*np.pi*f) / U
    K_e = 0.75 / L
    K_x_hat = K_x / K_e

    spl_h = 10 * np.log10(np.pow(M, 5) * (L*b)/(2*np.square(r_e)) \
                          * np.square(I) * np.square(rho_0) * np.pow(c_0, 4) \
                          * np.pow(K_x_hat, 3) / np.pow(1 + np.square(K_x_hat), 7/3)) + 78.4

    # Determine the low-frequency correction
    K_x_line = K_x * c / 2
    beta = np.sqrt(1 - np.square(M))

    S = np.sqrt(1 / (2 * np.pi * K_x_line / np.square(beta) + 1 / (1 + 2.4 * K_x_line / np.square(beta))))

    LFC = 10 * (np.square(S) * M * np.square(K_x_line) / np.square(beta)) * (1 + 9 * np.square(alpha))

    # Determine the corrected SPL
    spl_amiet = spl_h + 10 * np.log10(LFC / (1 + LFC))

    return spl_amiet

def moriarty(blade, f, r, c, U):
    """
    Calculate the SPL correction due to the airfoil thickness.
    
    Args:
        blade : Blade -- blade object
        f : np.array -- frequency array, [Hz]
        r : np.array -- radius array, [m]
        c : np.array -- chord array, [m]
        U : np.array -- velocity array, [m/s]

    Returns:
        delta_spl : np.array -- SPL correction array, [dB]
    """

    # Determine the number of frequencies and panels
    n_frequencies = f.size
    n_panels = U.size

    # Initialise the SPL matrix
    delta_spl = np.zeros((n_frequencies, n_panels))

    # Check for high Strouhal numbers
    # SPL will be low, so model accuracy is not critical
    St = (2*np.pi*f) * c / U

    if np.any(St > 75):
        print("High Strouhal numbers detected (St>75)!")

    # Iterate through each panel
    for i in range(n_panels):

        # Determine the neighbouring airfoils
        index_inboard  = blade.data[blade.data["pos"] < r[i]].index.max()
        index_outboard = blade.data[blade.data["pos"] > r[i]].index.min()

        r_inboard  = blade.data["pos"].iloc[index_inboard]
        r_outboard = blade.data["pos"].iloc[index_outboard]

        c_inboard  = blade.data["chord"].iloc[index_inboard]
        c_outboard = blade.data["chord"].iloc[index_outboard]

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

        airfoil_inboard  = None
        airfoil_outboard = None

        for airfoil in blade.airfoils:

            if airfoil.path == airfoil_path_inboard:
                airfoil_inboard = airfoil

            if airfoil.path == airfoil_path_outboard:
                airfoil_outboard = airfoil

        if airfoil_inboard is None or airfoil_outboard is None:
            print("Airfoil not found!")
            sys.exit(1)

        # Determine the thickness of the neighbouring airfoils
        xc_inboard  = airfoil_inboard.data["x/c"]
        yc_inboard  = airfoil_inboard.data["y/c"]
        xc_outboard = airfoil_outboard.data["x/c"]
        yc_outboard = airfoil_outboard.data["y/c"]

        xc_top_inboard  = xc_inboard[yc_inboard  > 0][::-1]
        yc_top_inboard  = yc_inboard[yc_inboard  > 0][::-1]
        xc_top_outboard = xc_outboard[yc_outboard > 0][::-1]
        yc_top_outboard = yc_outboard[yc_outboard > 0][::-1]

        xc_bottom_inboard  = xc_inboard[yc_inboard < 0]
        yc_bottom_inboard  = yc_inboard[yc_inboard < 0]
        xc_bottom_outboard = xc_outboard[yc_outboard < 0]
        yc_bottom_outboard = yc_outboard[yc_outboard < 0]

        yc_top_1_inboard     = np.interp(0.01, xc_top_inboard, yc_top_inboard)
        yc_top_1_outboard    = np.interp(0.01, xc_top_outboard, yc_top_outboard)
        yc_bottom_1_inboard  = np.interp(0.01, xc_bottom_inboard, yc_bottom_inboard)
        yc_bottom_1_outboard = np.interp(0.01, xc_bottom_outboard, yc_bottom_outboard)

        yc_top_10_inboard     = np.interp(0.10, xc_top_inboard, yc_top_inboard)
        yc_top_10_outboard    = np.interp(0.10, xc_top_outboard, yc_top_outboard)
        yc_bottom_10_inboard  = np.interp(0.10, xc_bottom_inboard, yc_bottom_inboard)
        yc_bottom_10_outboard = np.interp(0.10, xc_bottom_outboard, yc_bottom_outboard)

        t_1_inboard  = (yc_top_1_inboard - yc_bottom_1_inboard) * c_inboard
        t_1_outboard = (yc_top_1_outboard - yc_bottom_1_outboard) * c_outboard

        t_10_inboard  = (yc_top_10_inboard - yc_bottom_10_inboard) * c_inboard
        t_10_outboard = (yc_top_10_outboard - yc_bottom_10_outboard) * c_outboard

        # Determine the relative thickness of the interpolated airfoil
        D_rel_1  = np.interp(r[i], [r_inboard, r_outboard],  [t_1_inboard,  t_1_outboard]) / c[i]
        D_rel_10 = np.interp(r[i], [r_inboard, r_outboard], [t_10_inboard, t_10_outboard]) / c[i]

        # Determine the graph slope
        IT = D_rel_1 + D_rel_10
        SL = 1.123 * IT + 5.317 * np.square(IT)

        # Determine the SPL correction
        delta_spl[:, i] = -SL * (St[:, i] + 5)

    return delta_spl

def directivity(M, theta_e, phi_e):
    """
    Determine the directivity function.

    Arguments:
        M : np.array -- Mach number array, [-]
        theta_e : np.array -- angle array, [rad]
        phi_e : np.array -- angle array, [rad]

    Returns:
        D_L_line : np.array -- directivity array, [-]
    """

    # Determine the low-frequency directivity function 
    D_L_line = np.square(np.sin(theta_e)) * np.square(np.sin(phi_e)) \
               / np.pow(1 + M * np.cos(theta_e), 4)

    return D_L_line

def inflow_noise(f, turbine, results):
    """
    Determine the inflow noise.

    Arguments:
        f : np.array -- frequency array, [Hz]
        turbine : Turbine -- turbine object
        results : pd.DataFrame -- results dataframe

    Returns:
        spl : np.array -- SPL array, [dB]
    """

    # Determine the number of panels
    n_panels = turbine.attributes["NUMPANELS"]

    # TODO: WRITE COMMENT
    chords = np.array(turbine.blade.data["chord"])
    radiuses = np.array(turbine.blade.data["pos"])

    # Initialise inputs
    U     = np.zeros(n_panels)
    alpha = np.zeros(n_panels)
    r     = np.zeros(n_panels)

    for i in range(n_panels):
        U[i] = results[f"Total_Velocity_Blade_1_PAN_{i}_[m/s]"]
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
    r_e = 100* np.ones(n_panels) # TODO: Fix

    # Calculate SPL contributions
    spl_amiet = amiet(f, b, c, U, alpha, r_e, L, I, C_0, RHO_0)
    delta_spl = moriarty(turbine.blade, f, r, c, U)

    # Sum SPL contributions
    if SPL_CORRECTION:
        spl = spl_amiet + delta_spl + 10
    else:
        spl = spl_amiet + delta_spl

    return spl

# if __name__ == "__main__":

#     # Load data
#     turbine = Turbine("temp\\DTU_10MW\\DTU_10MW_RWT.trb")
#     results = pd.read_csv("temp\\results.txt", skiprows=2, delimiter="\t").iloc[TIMESTEP]

#     # Determine 1/3 band octave frequencies
#     f, _, _ = octave(F_MIN, F_MAX, F_REF, BASE_10)

#     f = f[:, np.newaxis]

#     # Determine the SPL
#     spl = inflow_noise(f, turbine, results)

#     # Plot spectra
#     for i in range(spl.shape[1]):
#         plt.plot(f, spl[:, i], label=f"Panel {i}")

#     plt.xlim(F_MIN, 1000)
#     plt.ylim(0, 80)
#     plt.xlabel("Frequency, [Hz]")
#     plt.ylabel("SPL, [dB]")
#     plt.xscale("log")
#     plt.yscale("linear")
#     plt.legend()
#     plt.grid(which="both")
#     plt.show()
