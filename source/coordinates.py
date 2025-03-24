"""
Author:   T. Moreira da Fonte Fonseca Teles
Email:    tmoreiradafont@tudelft.nl
Date:     2025-03-42
License:  GNU GPL 3.0

Coordinate system transformations.

Classes:
    None

Functions:
    translate
    rotate
    global_to_turbine
    turbine_to_global
    turbine_to_nacelle
    nacelle_to_turbine
    nacelle_to_hub
    hub_to_nacelle
    hub_to_blade
    blade_to_hub
    blade_to_airfoil
    airfoil_to_blade

Exceptions:
    None
"""

import sys

import numpy as np


def translate(p_x, p_y, p_z):
    """
    Returns the vector for a translation.

    Parameters:
        p_x : np.array -- translation in the x direction, [m]
        p_y : np.array -- translation in the y direction, [m]
        p_z : np.array -- translation in the z direction, [m]

    Returns:
        vector : np.array -- translation vector
    """

    if isinstance(p_x, np.ndarray):
        N = p_x.shape[0]

    elif isinstance(p_y, np.ndarray):
        N = p_y.shape

    elif isinstance(p_z, np.ndarray):
        N = p_z.shape

    else:
        N = 1

    p_x = np.resize(p_x, N)
    p_y = np.resize(p_y, N)
    p_z = np.resize(p_z, N)

    v_t = np.array([[p_x],
                    [p_y],
                    [p_z]])

    return np.transpose(v_t, (2, 0, 1))


def rotate(r_x, r_y, r_z, order):
    """
    Returns the matrix for a rotation.

    Parameters:
        r_x : np.array -- rotation angle about the x axis, [rad]
        r_y : np.array -- rotation angle about the y axis, [rad]
        r_z : np.array -- rotation angle about the z axis, [rad]
        order : str -- order of rotation

    Returns:
        matrix : np.array -- Nx3x3 rotation matrix
    """

    if isinstance(r_x, np.ndarray):
        N = r_x.shape

    elif isinstance(r_y, np.ndarray):
        N = r_y.shape

    elif isinstance(r_z, np.ndarray):
        N = r_z.shape

    else:
        N = 1

    r_x = np.resize(r_x, N)
    r_y = np.resize(r_y, N)
    r_z = np.resize(r_z, N)

    m_rot_x = np.array([[ np.ones(N), np.zeros(N),  np.zeros(N)],
                        [np.zeros(N), np.cos(r_x), -np.sin(r_x)],
                        [np.zeros(N), np.sin(r_x),  np.cos(r_x)]])

    m_rot_y = np.array([[ np.cos(r_y), np.zeros(N), np.sin(r_y)],
                        [ np.zeros(N),  np.ones(N), np.zeros(N)],
                        [-np.sin(r_y), np.zeros(N), np.cos(r_y)]])

    m_rot_z = np.array([[np.cos(r_z), -np.sin(r_z), np.zeros(N)],
                        [np.sin(r_z),  np.cos(r_z), np.zeros(N)],
                        [np.zeros(N),  np.zeros(N),  np.ones(N)]])

    m_rot_x = np.transpose(m_rot_x, (2, 0, 1))
    m_rot_y = np.transpose(m_rot_y, (2, 0, 1))
    m_rot_z = np.transpose(m_rot_z, (2, 0, 1))

    if order == "xyz":
        return m_rot_z @ m_rot_y @ m_rot_x

    elif order == "xzy":
        return m_rot_y @ m_rot_z @ m_rot_x

    elif order == "yxz":
        return m_rot_z @ m_rot_x @ m_rot_y

    elif order == "yzx":
        return m_rot_x @ m_rot_z @ m_rot_y

    elif order == "zxy":
        return m_rot_y @ m_rot_x @ m_rot_z

    elif order == "zyx":
        return m_rot_x @ m_rot_y @ m_rot_z

    else:
        print("Invalid order of rotation!")
        sys.exit(1)


def global_to_turbine(x_g, p_x, p_y, p_z, r_x, r_y, r_z):
    """
    Returns the turbine coordinates from the global coordinates.

    Parameters:
        x_g : np.array -- position in the global coordinate system, [m]
        p_x : float -- translation in the x direction, [m]
        p_y : float -- translation in the y direction, [m]
        p_z : float -- translation in the z direction, [m]
        r_x : float -- rotation angle about the x axis, [rad]
        r_y : float -- rotation angle about the y axis, [rad]
        r_z : float -- rotation angle about the z axis, [rad]

    Returns:
        x_t : np.array -- position in the turbine coordinate system, [m]
    """

    return np.transpose(rotate(r_x, r_y, r_z, "zyx"), (0, 2, 1)) @ (x_g - translate(p_x, p_y, p_z))


def turbine_to_global(x_t, p_x, p_y, p_z, r_x, r_y, r_z):
    """
    Returns the global coordinates from the turbine coordinates.

    Parameters:
        x_t : np.array -- position in the turbine coordinate system, [m]
        pos_x : float -- translation in the x direction, [m]
        pos_y : float -- translation in the y direction, [m]
        pos_z : float -- translation in the z direction, [m]
        rot_x : float -- rotation angle about the x axis, [rad]
        rot_y : float -- rotation angle about the y axis, [rad]
        rot_z : float -- rotation angle about the z axis, [rad]

    Returns:
        x_g : np.array -- position in the global coordinate system, [m]
    """

    return rotate(r_x, r_y, r_z, "zyx") @ x_t + translate(p_x, p_y, p_z)


def turbine_to_nacelle(x_t, tower_height, shaft_tilt, yaw):
    """
    Returns the nacelle coordinates from the turbine coordinates.

    Parameters:
        x_t : np.array -- position in the turbine coordinate system, [m]
        tower_height : float -- tower height, [m]
        shaft_tilt : float -- shaft tilt angle, [rad]
        yaw : float -- yaw angle, [rad]

    Returns:
        x_n : np.array -- position in the nacelle coordinate system, [m]
    """

    return np.transpose(rotate(0, shaft_tilt, yaw, "xyz"), (0, 2, 1)) @ (x_t - translate(0, 0, tower_height))


def nacelle_to_turbine(x_n, tower_height, shaft_tilt, yaw):
    """
    Returns the turbine coordinates from the nacelle coordinates.

    Parameters:
        x_n : np.array -- position in the nacelle coordinate system, [m]
        tower_height : float -- tower height, [m]
        shaft_tilt : float -- shaft tilt angle, [rad]
        yaw : float -- yaw angle, [rad]

    Returns:
        x_t : np.array -- position in the turbine coordinate system, [m]
    """

    return rotate(0, shaft_tilt, yaw, "xyz") @ x_n + translate(0, 0, tower_height)


def nacelle_to_hub(x_n, overhang, azimuth):
    """
    Returns the hub coordinates from the nacelle coordinates.

    Parameters:
        x_n : np.array -- position in the nacelle coordinate system, [m]
        overhang : float -- rotor overhang, [m]
        azimuth : float -- azimuthal rotor angle, [rad]

    Returns:
        x_h : np.array -- position in the hub coordinate system, [m]
    """

    return np.transpose(rotate(azimuth, 0, 0, "xyz"), (0, 2, 1)) \
           @ (x_n - translate(-overhang, 0, 0))


def hub_to_nacelle(x_h, overhang, azimuth):
    """
    Returns the nacelle coordinates from the hub coordinates.

    Parameters:
        x_h : np.array -- position in the hub coordinate system, [m]
        overhang : float -- rotor overhang, [m]
        azimuth : float -- azimuthal rotor angle, [rad]

    Returns:
        x_n : np.array -- position in the nacelle coordinate system, [m]
    """

    return rotate(azimuth, 0, 0, "xyz") @ x_h + translate(-overhang, 0, 0)


def hub_to_blade(x_h, cone, pitch):
    """
    Returns the blade coordinates from the hub coordinates.

    Parameters:
        x_h : np.array -- position in the hub coordinate system, [m]
        cone : float -- rotor cone angle, [rad]
        pitch : float -- collective blade pitch, [rad]
    """

    return np.transpose(rotate(0, -cone, -pitch, "xzy"), (0, 2, 1)) @ x_h


def blade_to_hub(x_b, cone, pitch):
    """
    Returns the hub coordinates from the blade coordinates.

    Parameters:
        x_b : np.array -- position in the blade coordinate system, [m]
        cone : float -- rotor cone angle, [rad]
        pitch : float -- collective blade pitch, [rad]

    Returns:
        x_h : np.array -- position in the hub coordinate system, [m]
    """

    return rotate(0, -cone, -pitch, "xzy") @ x_b


def blade_to_airfoil(x_b, pos, chord, twist, offset_x, offset_y, p_axis):
    """
    Returns the airfoil coordinates from the blade coordinates.

    Parameters:
        x_b : np.array -- position in the blade coordinate system, [m]
        pos : float -- radius of the blade, [m]
        chord : float -- chord of the blade, [m]
        twist : float -- twist angle of the blade, [rad]
        offset_x : float -- offset of the blade in the x-direction, [m]
        offset_y : float -- offset of the blade in the y-direction, [m]
        p_axis : float -- pitch axis, [-]

    Returns:
        x_a : np.array -- position in the airfoil coordinate system, [m]
    """

    return np.transpose(rotate(np.pi/2, 0, np.pi/2 - twist, "xyz"), (0, 2, 1)) \
           @ (x_b - translate(offset_y, offset_x, pos)) + translate(p_axis*chord, 0, 0)


def airfoil_to_blade(x_a, pos, chord, twist, offset_x, offset_y, p_axis):
    """
    Returns the blade coordinates from the airfoil coordinates.

    Parameters:
        x_a : np.array -- position in the airfoil coordinate system, [m]
        pos : float -- radius of the airfoil, [m]
        chord : float -- chord of the airfoil, [m]
        twist : float -- twist angle of the airfoil, [rad]
        offset_x : float -- offset of the airfoil in the x-direction, [m]
        offset_y : float -- offset of the airfoil in the y-direction, [m]
        p_axis : float -- pitch axis, [-]

    Returns:
        x_b : np.array -- position in the blade coordinate system, [m]
    """

    return rotate(np.pi/2, 0, np.pi/2 - twist, "xyz") @ (x_a - translate(p_axis*chord, 0, 0)) \
           + translate(offset_y, offset_x, pos)
