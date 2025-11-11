"""
Author:   T. Moreira da Fonte Fonseca Teles
Email:    tmoreiradafont@tudelft.nl
Date:     2025-11-11
License:  GNU GPL 3.0

Determine the displacement thickness.

Classes:
    None

Functions:
    displacement_thickness -- Determine the displacement thickness.

Exceptions:
    None
"""

import os

import numpy as np
import pandas as pd

from source.XFOIL.xfoil import XFoil


def displacement_thickness(blade, c, r, U, alpha, nu, xfoil_path, max_iter, \
                           x_c_upper, x_c_lower, n_crit, probe_upper, probe_lower):
    """
    Determine the displacement thickness.

    Arguments:
        blade : Blade -- blade object
        c : np.ndarray -- chord, [m]
        r : np.ndarray -- radius, [m]
        U : np.ndarray -- velocity, [m/s]
        alpha : np.ndarray -- angle of attack, [rad]
        nu : float -- kinematic viscosity, [m^2/s]
        xfoil_path : str -- path to the XFOIL executable
        max_iter : int -- maximum number of XFOIL iterations, [-]
        x_c_upper : float -- transition position on the upper surface, [-]
        x_c_lower : float -- transition position on the lower surface, [-]
        n_crit : float -- critical amplification factor, [-]
        probe_upper : float -- probe position on the upper surface, [-]
        probe_lower : float -- probe position on the lower surface, [-]

    Returns:
        delta_star_upper : np.ndarray -- displacement thickness on the upper surface, [m]
        delta_star_lower : np.ndarray -- displacement thickness on the upper surface, [m]
    """

    # Determine the Reynolds number
    re = U * c / nu

    # Initialise the displacement thickness arrays
    delta_star_upper = np.empty(r.shape)
    delta_star_lower = np.empty(r.shape)

    # Iterate over the blade panels
    for i in range(len(r)):

        # Determine the airfoil indices
        index_0 = np.searchsorted(blade.radius, r[i]) - 1
        index_1 = np.searchsorted(blade.radius, r[i])

        # Determine the airfoil radiuses
        radius_0 = blade.radius[index_0]
        radius_1 = blade.radius[index_1]

        # Determine the airfoil paths
        path_0 = blade.polar[index_0].airfoil.path
        path_1 = blade.polar[index_1].airfoil.path

        # Determine the interpolation fraction
        fraction = (r[i] - radius_0) / (radius_1 - radius_0)

        # Determine the current working directory
        cwd = os.path.dirname(path_0)

        # Determine the output path
        path_out = os.path.join(cwd, f"xfoil_{i:02d}.csv")

        # Determine the relative paths
        name_0 = os.path.basename(path_0)
        name_1 = os.path.basename(path_1)
        name_out = os.path.basename(path_out)

        # Initialise XFOIL
        xfoil = XFoil(xfoil_path, cwd)

        # Run the XFOIL simulation
        xfoil.run(name_0, name_1, name_out, fraction, re[i], alpha[i], max_iter, x_c_upper, x_c_lower, n_crit)

        # Read the output file
        data = pd.read_csv(path_out, delimiter=r"\s+", names=["s/c", "x/c", "y/c", \
                           "U_e/U", "delta_star/c", "theta/c", "C_f", "H"], skiprows=1)

        # Filter and sort the boundary layer data
        index = data["x/c"].idxmin()
        data_upper = data[data["x/c"] <= 1.0].iloc[:index + 1][::-1]
        data_lower = data[data["x/c"] <= 1.0].iloc[index:]

        # Determine the displacement thickness at the probe locations
        delta_star_upper[i] = np.interp(probe_upper, data_upper["x/c"], data_upper["delta_star/c"]) * c[i]
        delta_star_lower[i] = np.interp(probe_lower, data_lower["x/c"], data_lower["delta_star/c"]) * c[i]

        # Remove the output file
        os.remove(path_out)

    return delta_star_upper, delta_star_lower
