"""
Author:   T. Moreira da Fonte Fonseca Teles
Email:    tmoreiradafont@tudelft.nl
Date:     2025-09-17
License:  GNU GPL 3.0

Determine the coordinate system transformations.

Classes:
    None

Functions:
    transform
    nacelle_to_turbine
    hub_to_nacelle
    blade_to_hub
    airfoil_to_blade
    freestream_to_airfoil

Exceptions:
    None
"""

import sys

import numpy as np

from source.settings import SIMULATION_PATH
from source.QBlade.simulation import Simulation


def transform(t_x, t_y, t_z, r_x, r_y, r_z, order):
    """
    Determine the transformation matrix.

    Parameters:
        t_x : float -- translation in the x direction, [m]
        t_y : float -- translation in the y direction, [m]
        t_z : float -- translation in the z direction, [m]
        r_x : float -- rotation about the x-axis, [rad]
        r_y : float -- rotation about the y-axis, [rad]
        r_z : float -- rotation about the z-axis, [rad]
        order : str -- rotation order
    """

    # Determine the rotation matrix about the x-axis
    rotation_x = np.array([[1.0,         0.0,          0.0, 0.0],
                           [0.0, np.cos(r_x), -np.sin(r_x), 0.0],
                           [0.0, np.sin(r_x),  np.cos(r_x), 0.0],
                           [0.0,         0.0,          0.0, 1.0]])

    # Determine the rotation matrix about the y-axis
    rotation_y = np.array([[ np.cos(r_y), 0.0, np.sin(r_y), 0.0],
                           [         0.0, 1.0,         0.0, 0.0],
                           [-np.sin(r_y), 0.0, np.cos(r_y), 0.0],
                           [         0.0, 0.0,         0.0, 1.0]])

    # Determine the rotation matrix about the z axis
    rotation_z = np.array([[np.cos(r_z), -np.sin(r_z), 0.0, 0.0],
                           [np.sin(r_z),  np.cos(r_z), 0.0, 0.0],
                           [        0.0,          0.0, 1.0, 0.0],
                           [        0.0,          0.0, 0.0, 1.0]])

    # Determine the translation matrix
    translation = np.array([[1.0, 0.0, 0.0, t_x],
                            [0.0, 1.0, 0.0, t_y],
                            [0.0, 0.0, 1.0, t_z],
                            [0.0, 0.0, 0.0, 1.0]])

    # Determine the transformation matrix
    rotation = {
        "xyz": rotation_z @ rotation_y @ rotation_x,
        "xzy": rotation_y @ rotation_z @ rotation_x,
        "yxz": rotation_z @ rotation_x @ rotation_y,
        "yzx": rotation_x @ rotation_z @ rotation_y,
        "zxy": rotation_y @ rotation_x @ rotation_z,
        "zyx": rotation_x @ rotation_y @ rotation_z,
        }

    if order not in rotation:
        print("Invalid order of rotation!")
        sys.exit(1)

    return translation @ rotation[order]


def nacelle_to_turbine(tower_height, shaft_tilt, yaw):
    """
    Determine the transformation matrix from the nacelle
    coordinate system to the turbine coordinate system.

    Parameters:
        tower_height : float -- tower height, [m]
        shaft_tilt : float -- shaft tilt angle, [rad]
        yaw : float -- yaw angle, [rad]

    Returns:
        matrix : np.ndarray -- transformation matrix
    """

    return transform(0.0, 0.0, tower_height, 0.0, shaft_tilt, yaw, "xyz")


def hub_to_nacelle(rotor_overhang, azimuth):
    """
    Determine the transformation matrix from the hub
    coordinate system to the nacelle coordinate system.

    Parameters:
        rotor_overhang : float -- rotor overhang, [m]
        azimuth : float -- azimuth angle, [rad]

    Returns:
        matrix : np.ndarray -- transformation matrix
    """

    return transform(-rotor_overhang, 0.0, 0.0, azimuth, 0.0, 0.0, "xyz")


def blade_to_hub(rotor_cone, pitch):
    """
    Determine the transformation matrix from the blade
    coordinate system to the hub coordinate system.

    Parameters:
        rotor_cone : float -- rotor cone angle, [rad]
        pitch : float -- pitch angle, [rad]

    Returns:
        matrix : np.ndarray -- transformation matrix
    """

    return transform(0.0, 0.0, 0.0, 0.0, -rotor_cone, -pitch, "xzy")


def airfoil_to_blade(radius, chord, twist, offset_x, offset_y, pitch_axis, x_c):
    """
    Determine the transformation matrix from the airfoil
    coordinate system to the blade coordinate system.

    Parameters:
        radius : float -- spanwise position, [m]
        chord : float -- chord, [m]
        twist : float -- twist angle, [rad]
        offset_x : float -- offset in the x-direction, [m]
        offset_y : float -- offset in the y-direction, [m]
        pitch_axis : float -- pitch axis position, [-]
        x_c : float -- origin position, [-]

    Returns:
        matrix : np.ndarray -- transformation matrix
    """

    return transform(offset_y, offset_x, radius, np.pi/2.0, 0.0, np.pi/2.0 - twist, "xyz") \
           @ transform((x_c - pitch_axis) * chord, 0.0, 0.0, 0.0, 0.0, 0.0, "xyz")


def freestream_to_airfoil(alpha):
    """
    Determine the transformation matrix from the freestream
    coordinate system to the airfoil coordinate system.

    Parameters:
        alpha : float -- angle of attack, [rad]

    Returns:
        matrix : np.ndarray -- transformation matrix
    """

    return transform(0.0, 0.0, 0.0, 0.0, -alpha, 0.0, "xyz")

if __name__ == "__main__":

    # Open the output file
    f = open("coordinates.obj", "w", encoding="utf-8")

    # Create the Simulation object
    simulation = Simulation(SIMULATION_PATH)
    turbine = simulation.turbine
    blade = turbine.blade

    # Iterate through all blades
    for i in range(turbine.n_blades):

        # Determine the azimuth angle
        azimuth = 2.0 * np.pi * i / turbine.n_blades

        # Iterate through all panels
        for j in range(len(blade.radius)):

            # Determine the blade properties
            radius = blade.radius[j]
            chord = blade.chord[j]
            twist = blade.twist[j]
            offset_x = blade.offset_x[j]
            offset_y = blade.offset_y[j]
            pitch_axis = blade.pitch_axis[j]
            airfoil = blade.polar[j].airfoil

            # Determine the airfoil coordinates
            x = airfoil.coordinates["x/c"] * chord
            y = np.zeros(len(x))
            z = airfoil.coordinates["y/c"] * chord
            w = np.ones(len(x))

            # Determine the transformation matrices
            matrix_ab = airfoil_to_blade(radius, chord, twist, offset_x, offset_y, pitch_axis, 0.0)
            matrix_bh = blade_to_hub(turbine.rotor_cone, 0.0)
            matrix_hn = hub_to_nacelle(turbine.rotor_overhang, azimuth)
            matrix_nt = nacelle_to_turbine(turbine.tower_height, turbine.shaft_tilt, 0.0)

            # Transform the airfoil coordinates
            x_t = matrix_nt @ matrix_hn @ matrix_bh @ matrix_ab @ np.array([x, y, z, w])

            # Write the airfoil coordinates
            for k in range(x_t.shape[1]):
                f.write(f"v {x_t[0, k]} {x_t[1, k]} {x_t[2, k]}\n")
