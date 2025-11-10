"""
Author:   T. Moreira da Fonte Fonseca Teles
Email:    tmoreiradafont@tudelft.nl
Date:     2025-11-10
License:  GNU GPL 3.0

Determine the boundary layer displacement thickness using XFOIL.

Classes:
    None

Functions:
    displacement_thickness -- Determine the boundary layer displacement thickness using XFOIL.

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
    Determine the boundary layer displacement thickness using XFOIL.

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
        probe_upper : np.ndarray -- probe position on the upper surface, [-]
        probe_lower : np.ndarray -- probe position on the lower surface, [-]

    Returns:
        delta_star_upper : np.ndarray -- upper displacement thickness, [m]
        delta_star_lower : np.ndarray -- lower displacement thickness, [m]
    """

    # Determine the Reynolds number
    re = U * c / nu

    # Determine the airfoil indices
    index_0 = np.searchsorted(blade.radius, r) - 1
    index_1 = np.searchsorted(blade.radius, r)

    # Initialise the displacement thickness arrays
    delta_star_upper = np.empty(r.shape)
    delta_star_lower = np.empty(r.shape)

    # Iterate over the blade panels
    for i in range(len(r)):

        # Determine the airfoil paths
        path_0 = blade.polar[index_0].airfoil.path
        path_1 = blade.polar[index_1].airfoil.path

        # Determine the current working directory
        cwd = os.path.dirname(path_0)

        # Determine the relative airfoil paths
        path_rel_0 = os.path.basename(path_0)
        path_rel_1 = os.path.basename(path_1)

        # Determine the output dump paths
        path_rel_out = f"xfoil_dump_{i:2d}.csv"
        path_out = os.path.join(cwd, path_out)

        # Determine the interpolation fraction
        radius_0 = blade.radius[index_0]
        radius_1 = blade.radius[index_1]

        fraction = (r[i] - radius_0) / (radius_1 - radius_0)

        # Initialise XFOIL
        xfoil = XFoil(xfoil_path, cwd)

        # Run the XFOIL simulation
        xfoil.run(path_0, path_1, fraction, path_rel_out, re[i], alpha[i], max_iter, x_c_upper, x_c_lower, n_crit)

        # Read the output file
        data = pd.read_csv(path_out, delimiter=r"\s+", names=["s/c", "x/c", "y/c", "U_e/U", \
                           "delta_star/c", "theta/c", "C_f", "H"], skiprows=1)

        # Filter and sort the boundary layer data
        index = data["x/c"].idxmin()
        upper = data[data["x/c"] <= 1.0].iloc[:index + 1][::-1]
        lower = data[data["x/c"] <= 1.0].iloc[index:]

        # Determine the displacement thickness at the probe locations
        delta_star_upper[i] = np.interp(probe_upper, upper["x/c"], upper["delta_star/c"]) * c[i]
        delta_star_lower[i] = np.interp(probe_lower, lower["x/c"], lower["delta_star/c"]) * c[i]

        # Remove the output file
        os.remove(path_out)

    return delta_star_upper, delta_star_lower
