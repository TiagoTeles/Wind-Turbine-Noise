"""
Author:   T. Moreira da Fonte Fonseca Teles
Email:    tmoreiradafont@tudelft.nl
Date:     2025-03-05
License:  GNU GPL 3.0

Miscellaneous functions.

Classes:
    None

Functions:
    octave
    E
    translate
    rotate

Exceptions:
    None
"""

import numpy as np
import scipy as sp


def octave(f_min, f_max, f_ref, base_10=True):
    """
    Determine the center, lower, and upper frequencies of the 1/3 octave band.

    Arguments:
        f_min : float -- Minimum frequency, [Hz]
        f_max : float -- Maximum frequency, [Hz]
        f_ref : float -- Reference frequency, [Hz]
        base_10 : bool -- Use base 10?

    Returns:
        f_center : np.array -- Center frequencies, [Hz]
        f_lower : np.array -- Lower frequencies, [Hz]
        f_upper : np.array -- Upper frequencies, [Hz]
    """

    if base_10:
        # Determine smallest and largest index
        min_index = np.floor(10 * np.log10(f_min / f_ref) + 0.5)
        max_index =  np.ceil(10 * np.log10(f_max / f_ref) - 0.5)

        # Determine center, lower and upper frequencies
        f_center = f_ref * np.pow(10, np.arange(min_index, max_index+1) / 10)
        f_lower = f_center / np.pow(10, 1/20)
        f_upper = f_center * np.pow(10, 1/20)

    else:
        # Determine smallest and largest index
        min_index = np.floor(3 * np.log2(f_min / f_ref) + 0.5)
        max_index =  np.ceil(3 * np.log2(f_max / f_ref) - 0.5)

        # Determine center, lower and upper frequencies
        f_center = f_ref * np.pow(2, np.arange(min_index, max_index+1) / 3)
        f_lower = f_center / np.pow(2, 1/6)
        f_upper = f_center * np.pow(2, 1/6)

    return f_center, f_lower, f_upper

def E(x):
    """
    Determine the combination of Fresnel integrals defined in Roger and Moreau (2005).

    Arguments:
        x : np.array -- argument, [-]

    Returns:
        E : np.array -- C_2 - i * S_2, [-]
    """

    s_2, c_2 = sp.special.fresnel(np.sqrt(2*x/np.pi))

    return c_2 - 1j * s_2

def translate(pos_x, pos_y, pos_z):
    """
    Returns the vector for a translation.

    Parameters:
        pos_x : float -- translation in the x direction, [m]
        pos_y : float -- translation in the y direction, [m]
        pos_z : float -- translation in the z direction, [m]
    """

    return np.array([pos_x, pos_y, pos_z])[:, np.newaxis]

def rotate(rot_x, rot_y, rot_z, order):
    """
    Returns the transformation matrix for a rotation.

    Parameters:
        rot_x : float -- rotation angle about the x axis, [rad]
        rot_y : float -- rotation angle about the y axis, [rad]
        rot_z : float -- rotation angle about the z axis, [rad]
        order : str -- order of rotation

    Returns:
        matrix : np.array -- 3x3 rotation matrix
    """

    m_rot_x = np.array([[1,             0,              0],
                        [0, np.cos(rot_x), -np.sin(rot_x)],
                        [0, np.sin(rot_x),  np.cos(rot_x)]])

    m_rot_y = np.array([[ np.cos(rot_y), 0, np.sin(rot_y)],
                        [            0,  1,             0],
                        [-np.sin(rot_y), 0, np.cos(rot_y)]])

    m_rot_z = np.array([[np.cos(rot_z), -np.sin(rot_z), 0],
                        [np.sin(rot_z),  np.cos(rot_z), 0],
                        [            0,              0, 1]])

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
