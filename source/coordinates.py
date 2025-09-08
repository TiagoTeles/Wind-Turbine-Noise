"""
Author:   T. Moreira da Fonte Fonseca Teles
Email:    tmoreiradafont@tudelft.nl
Date:     2025-09-08
License:  GNU GPL 3.0

Determine the coordinate system transformations.

Classes:
    None

Functions:
    transform
    turbine_to_global
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
    if order == "xyz":
        return translation @ rotation_z @ rotation_y @ rotation_x

    elif order == "xzy":
        return translation @ rotation_y @ rotation_z @ rotation_x

    elif order == "yxz":
        return translation @ rotation_z @ rotation_x @ rotation_y

    elif order == "yzx":
        return translation @ rotation_x @ rotation_z @ rotation_y

    elif order == "zxy":
        return translation @ rotation_y @ rotation_x @ rotation_z

    elif order == "zyx":
        return translation @ rotation_x @ rotation_y @ rotation_z

    else:
        print("Invalid order of rotation!")
        sys.exit(1)


def turbine_to_global(t_x, t_y, t_z, r_x, r_y, r_z):
    """
    Determine the transformation matrix from the turbine
    coordinate system to the global coordinate system.

    Parameters:
        t_x : float -- translation in the x direction, [m]
        t_y : float -- translation in the y direction, [m]
        t_z : float -- translation in the z direction, [m]
        r_x : float -- rotation about the x-axis, [rad]
        r_y : float -- rotation about the y-axis, [rad]
        r_z : float -- rotation about the z-axis, [rad]

    Returns:
        matrix : np.array -- transformation matrix
    """

    return transform(t_x, t_y, t_z, r_x, r_y, r_z, "zyx")


def nacelle_to_turbine(tower_height, shaft_tilt_angle, phi):
    """
    Determine the transformation matrix from the nacelle
    coordinate system to the turbine coordinate system.

    Parameters:
        tower_height : float -- tower height, [m]
        shaft_tilt_angle : float -- shaft tilt angle, [rad]
        phi : float -- yaw angle, [rad]

    Returns:
        matrix : np.array -- transformation matrix
    """

    return transform(0.0, 0.0, tower_height, 0.0, shaft_tilt_angle, phi, "xyz")


def hub_to_nacelle(rotor_overhang, psi):
    """
    Determine the transformation matrix from the hub
    coordinate system to the nacelle coordinate system.

    Parameters:
        rotor_overhang : float -- rotor overhang, [m]
        psi : float -- azimuth angle, [rad]

    Returns:
        matrix : np.array -- transformation matrix
    """

    return transform(-rotor_overhang, 0.0, 0.0, psi, 0.0, 0.0, "xyz")


def blade_to_hub(rotor_cone_angle, theta):
    """
    Determine the transformation matrix from the blade
    coordinate system to the hub coordinate system.

    Parameters:
        rotor_cone_angle : float -- rotor cone angle, [rad]
        theta : float -- blade pitch angle, [rad]

    Returns:
        matrix : np.array -- transformation matrix
    """

    return transform(0.0, 0.0, 0.0, 0.0, -rotor_cone_angle, -theta, "xzy")


def airfoil_to_blade(radius, chord, twist_angle, offset_x, offset_y, pitch_axis, x_c):
    """
    Determine the transformation matrix from the airfoil
    coordinate system to the blade coordinate system.

    Parameters:
        radius : float -- spanwise position, [m]
        chord : float -- chord, [m]
        twist_angle : float -- twist angle, [rad]
        offset_x : float -- offset in the x-direction, [m]
        offset_y : float -- offset in the y-direction, [m]
        pitch_axis : float -- pitch axis position, [-]
        x_c : float -- origin position, [-]

    Returns:
        matrix : np.array -- transformation matrix
    """

    return transform(offset_y, offset_x, radius, np.pi/2.0, 0.0, np.pi/2.0 - twist_angle, "xyz") \
           @ transform((x_c - pitch_axis) * chord, 0.0, 0.0, 0.0, 0.0, 0.0, "xyz")


def freestream_to_airfoil(alpha):
    """
    Determine the transformation matrix from the freestream
    coordinate system to the airfoil coordinate system.

    Parameters:
        alpha : float -- angle of attack, [rad]

    Returns:
        matrix : np.array -- transformation matrix
    """

    return transform(0.0, 0.0, 0.0, 0.0, -alpha, 0.0, "xyz")
