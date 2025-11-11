"""
Author:   T. Moreira da Fonte Fonseca Teles
Email:    tmoreiradafont@tudelft.nl
Date:     2025-11-11
License:  GNU GPL 3.0

Determine the retarded source and observer coordinates.

Classes:
    None

Functions:
    retarded_coordinates

Exceptions:
    None
"""

import numpy as np

from source.coordinates.transform import nacelle_to_turbine, hub_to_nacelle, blade_to_hub
from source.coordinates.transform import airfoil_to_blade, freestream_to_airfoil


def retarded_coordinates(simulation, x_o_t, azimuth, radius, c_0):
    """
    Determine the retarded source and observer coordinates.

    Parameters:
        simulation : Simulation -- QBlade simulation
        x_o_t : np.ndarray -- observer coordinates in the turbine coordinate system, [m]
        azimuth : float -- azimuth angle at reception, [rad]
        radius : np.ndarray -- radius, [m]
        c_0 : float -- speed of sound, [m/s]

    Returns:
        x_s_t : np.ndarray -- source coordinates in the turbine coordinate system, [m]
        x_o_f : np.ndarray -- observer coordinates in the freestream coordinate system, [m]
    """

    # Determine the turbine and blade
    turbine = simulation.turbine
    blade = turbine.blade

    # Determine the number of sources and observers
    n_sources = radius.shape[0]
    n_observers = x_o_t.shape[1]

    # Determine the turbine properties
    tower_height = turbine.tower_height
    shaft_tilt = turbine.shaft_tilt
    rotor_overhang = turbine.rotor_overhang
    rotor_cone = turbine.rotor_cone

    # Determine the blade properties
    chord = blade.get_geometry("chord", radius)
    twist = blade.get_geometry("twist", radius)
    offset_x = blade.get_geometry("offset_x", radius)
    offset_y = blade.get_geometry("offset_y", radius)
    pitch_axis = blade.get_geometry("pitch_axis", radius)

    # Initialise the iteration parameters
    azimuth_e = np.ones((n_sources, n_observers)) * azimuth
    error = np.ones((n_sources, n_observers))

    # Initialise the coordinate arrays
    x_s_t = np.empty((3, n_sources, n_observers))
    x_o_f = np.empty((3, n_sources, n_observers))

    # Reshape the observer coordinates
    x_o_t = x_o_t[:, np.newaxis, :]

    # Iterate until convergence
    while np.any(np.abs(error) > 1.0E-3):

        # Determine the turbine operating conditions
        yaw = turbine.get_results("yaw", azimuth_e)
        pitch = turbine.get_results("pitch", azimuth_e)

        # Determine the source coordinates
        for i in range(n_sources):
            for j in range(n_observers):

                # Determine the blade operating conditions
                angle_of_attack = blade.get_results("angle_of_attack", azimuth_e[i, j], radius[i])

                # Determine the invividual transformation matrices
                matrix_nt = nacelle_to_turbine(tower_height, shaft_tilt, yaw[i, j])
                matrix_hn = hub_to_nacelle(rotor_overhang, azimuth_e[i, j])
                matrix_bh = blade_to_hub(rotor_cone, pitch[i, j])
                matrix_ab = airfoil_to_blade(radius[i], chord[i], twist[i], offset_x[i], offset_y[i], pitch_axis[i], 0.5)
                matrix_fa = freestream_to_airfoil(angle_of_attack)

                # Determine the complete transformation matrix
                matrix_ft = matrix_nt @ matrix_hn @ matrix_bh @ matrix_ab @ matrix_fa

                # Determine the source coordinates
                x_source = matrix_ft @ np.array([0, 0, 0, 1])

                # Remove the w coordinate
                x_s_t[:, i, j] = x_source[0:3]

        # Determine the source-observer distance
        distance = np.sqrt(np.sum(np.square(x_o_t - x_s_t), axis=0))

        # Determine the rotor angular velocity
        inflow_velocity = simulation.get_results("inflow_velocity", azimuth_e)
        tip_speed_ratio = turbine.get_results("tip_speed_ratio", azimuth_e)
        angular_velocity = tip_speed_ratio * inflow_velocity / blade.radius[-1]

        # Check if the solution has converged
        error = (distance / c_0) - ((azimuth - azimuth_e) / angular_velocity)

        # Determine the azimuth angle at emission
        azimuth_e = azimuth - angular_velocity * (distance / c_0)

    # Determine the observer coordinates
    for i in range(n_sources):
        for j in range(n_observers):

            # Determine the blade operating conditions
            angle_of_attack = blade.get_results("angle_of_attack", azimuth_e[i, j], radius[i])

            # Determine the invividual transformation matrices
            matrix_nt = nacelle_to_turbine(tower_height, shaft_tilt, yaw[i, j])
            matrix_hn = hub_to_nacelle(rotor_overhang, azimuth_e[i, j])
            matrix_bh = blade_to_hub(rotor_cone, pitch[i, j])
            matrix_ab = airfoil_to_blade(radius[i], chord[i], twist[i], offset_x[i], offset_y[i], pitch_axis[i], 0.5)
            matrix_fa = freestream_to_airfoil(angle_of_attack)

            # Determine the complete transformation matrix
            matrix_ft = matrix_nt @ matrix_hn @ matrix_bh @ matrix_ab @ matrix_fa

            # Determine the observer coordinates
            x_observer = np.linalg.inv(matrix_ft) @ np.append(x_o_t[:, :, j], 1.0)

            # Remove the w coordinate
            x_o_f[:, i, j] = x_observer[0:3]

    return x_s_t, x_o_f
