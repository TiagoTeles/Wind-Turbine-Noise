"""
Author:   T. Moreira da Fonte Fonseca Teles
Email:    tmoreiradafont@tudelft.nl
Date:     2025-05-12
License:  GNU GPL 3.0

Coordinate system transformations.

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


def transform(p_x, p_y, p_z, r_x, r_y, r_z, order):
    """
    Returns the transformation matrix.

    Parameters:
        p_x : float -- translation in the x direction, [m]
        p_y : float -- translation in the y direction, [m]
        p_z : float -- translation in the z direction, [m]
        r_x : float -- rotation about the x axis, [rad]
        r_y : float -- rotation about the y axis, [rad]
        r_z : float -- rotation about the z axis, [rad]
        order : str -- order of the rotation
    """

    # Determine the x rotation matrix
    m_rot_x = np.array([[1,           0,            0, 0],
                        [0, np.cos(r_x), -np.sin(r_x), 0],
                        [0, np.sin(r_x),  np.cos(r_x), 0],
                        [0,            0,           0, 1]])

    # Determine the y rotation matrix
    m_rot_y = np.array([[ np.cos(r_y), 0, np.sin(r_y), 0],
                        [           0, 1,           0, 0],
                        [-np.sin(r_y), 0, np.cos(r_y), 0],
                        [           0, 0,           0, 1]])

    # Determine the z rotation matrix
    m_rot_z = np.array([[np.cos(r_z), -np.sin(r_z), 0, 0],
                        [np.sin(r_z),  np.cos(r_z), 0, 0],
                        [          0,            0, 1, 0],
                        [          0,            0, 0, 1]])

    # Determine the translation matrix
    m_trans = np.array([[1, 0, 0, p_x],
                        [0, 1, 0, p_y],
                        [0, 0, 1, p_z],
                        [0, 0, 0,   1]])

    # Determine the transformation matrix
    if order == "xyz":
        return m_trans @ m_rot_z @ m_rot_y @ m_rot_x

    elif order == "xzy":
        return m_trans @ m_rot_y @ m_rot_z @ m_rot_x

    elif order == "yxz":
        return m_trans @ m_rot_z @ m_rot_x @ m_rot_y

    elif order == "yzx":
        return m_trans @ m_rot_x @ m_rot_z @ m_rot_y

    elif order == "zxy":
        return m_trans @ m_rot_y @ m_rot_x @ m_rot_z

    elif order == "zyx":
        return m_trans @ m_rot_x @ m_rot_y @ m_rot_z

    else:
        print("Invalid order of rotation!")
        sys.exit(1)


def turbine_to_global(p_x, p_y, p_z, r_x, r_y, r_z):
    """
    Returns the transformation matrix from the turbine 
    coordinate system to the global coordinate system.

    Parameters:
        p_x : float -- translation in the x direction, [m]
        p_y : float -- translation in the y direction, [m]
        p_z : float -- translation in the z direction, [m]
        r_x : float -- rotation about the x axis, [rad]
        r_y : float -- rotation about the y axis, [rad]
        r_z : float -- rotation about the z axis, [rad]

    Returns:
        matrix : np.array -- transformation matrix
    """

    return transform(p_x, p_y, p_z, r_x, r_y, r_z, "zyx")


def nacelle_to_turbine(tower_height, shaft_tilt, yaw):
    """
    Returns the transformation matrix from the nacelle 
    coordinate system to the turbine coordinate system.

    Parameters:
        tower_height : float -- tower height, [m]
        shaft_tilt : float -- shaft tilt angle, [rad]
        yaw : float -- yaw angle, [rad]

    Returns:
        matrix : np.array -- transformation matrix
    """

    return transform(0, 0, tower_height, 0, shaft_tilt, yaw, "xyz")


def hub_to_nacelle(overhang, azimuth):
    """
    Returns the transformation matrix from the hub 
    coordinate system to the nacelle coordinate system.

    Parameters:
        overhang : float -- rotor overhang, [m]
        azimuth : float -- azimuthal rotor angle, [rad]

    Returns:
        matrix : np.array -- transformation matrix
    """

    return transform(-overhang, 0, 0, azimuth, 0, 0, "xyz")


def blade_to_hub(cone, pitch):
    """
    Returns the transformation matrix from the blade 
    coordinate system to the hub coordinate system.

    Parameters:
        cone : float -- rotor cone angle, [rad]
        pitch : float -- collective blade pitch, [rad]

    Returns:
        matrix : np.array -- transformation matrix
    """

    return transform(0, 0, 0, 0, -cone, -pitch, "xzy")


def airfoil_to_blade(pos, chord, twist, offset_x, offset_y, p_axis, x_c):
    """
    Returns the transformation matrix from the airfoil 
    coordinate system to the blade coordinate system.

    Parameters:
        pos : float -- radius of the airfoil, [m]
        chord : float -- chord of the airfoil, [m]
        twist : float -- twist angle of the airfoil, [rad]
        offset_x : float -- offset of the airfoil in the x-direction, [m]
        offset_y : float -- offset of the airfoil in the y-direction, [m]
        p_axis : float -- pitch axis, [-]
        x_c : float -- coordinate system orgin, [-]

    Returns:
        matrix : np.array -- transformation matrix
    """

    return transform(offset_y, offset_x, pos, np.pi/2, 0, np.pi/2 - twist, "xyz") \
           @ transform((x_c-p_axis) * chord, 0, 0, 0, 0, 0, "xyz")


def freestream_to_airfoil(alpha):
    """
    Returns the transformation matrix from the freestream 
    coordinate system to the airfoil coordinate system.

    Parameters:
        alpha : float -- angle of attack, [rad]

    Returns:
        matrix : np.array -- transformation matrix
    """

    return transform(0, 0, 0, 0, -alpha, 0, "xyz")
