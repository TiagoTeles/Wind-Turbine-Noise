""" 
Author:   T. Moreira da Fonte Fonseca Teles
Email:    tmoreiradafont@tudelft.nl
Date:     2025-05-30
License:  GNU GPL 3.0

Read the data from .shd files.

Classes:
    None

Functions:
    read_shd

Exceptions:
    None
"""

import numpy as np


def read_shd(filename):
    """
    Read the .shd file.

    Arguments:
        filename : str -- path to the .shd file

    Returns:
        title : str -- title of the file
        frequency : np.ndarray -- frequency, [Hz]
        pressure : np.ndarray -- pressure, [Pa]
        s_x : np.ndarray -- source x coordinates, [m]
        s_y : np.ndarray -- source y coordinates, [m]
        s_z : np.ndarray -- source z coordinates, [m]
        r_z : np.ndarray -- receiver z coordinates, [m]
        r_r : np.ndarray -- receiver ranges, [m]
        r_theta : np.ndarray -- receiver bearings, [rad]
    """

    # Open the file
    f = open(filename, "rb")

    # Read the first record
    record_length = 4 * np.fromfile(f, np.int32, 1)[0]
    title = f.read(80)

    # Read the second record
    f.seek(1 * record_length)
    plot_type = f.read(10)

    # Read the third record
    f.seek(2 * record_length)
    n_frequencies = np.fromfile(f, np.int32, 1)[0]
    n_r_theta = np.fromfile(f, np.int32, 1)[0]
    n_s_x = np.fromfile(f, np.int32, 1)[0]
    n_s_y = np.fromfile(f, np.int32, 1)[0]
    n_s_z = np.fromfile(f, np.int32, 1)[0]
    n_r_z = np.fromfile(f, np.int32, 1)[0]
    n_r_r = np.fromfile(f, np.int32, 1)[0]
    frequency_0 = np.fromfile(f, np.float64, 1)[0]
    attenuation = np.fromfile(f, np.float64, 1)[0]

    # Read the fourth record
    f.seek(3 * record_length)
    frequency = np.fromfile(f, np.float64, n_frequencies)

    # Read the fifth record
    f.seek(4 * record_length)
    r_theta = np.fromfile(f, np.float64, n_r_theta)

    if plot_type[0 : 1] != "TL":

        # Read the sixth record
        f.seek(5 * record_length)
        s_x = np.fromfile(f, np.float64, n_s_x)

        # Read the seventh record
        f.seek(6 * record_length)
        s_y = np.fromfile(f, np.float64, n_s_y)

    else:

        # Read the sixth record
        f.seek(5 * record_length)
        s_x = np.fromfile(f, np.float64, 2)
        s_x = np.linspace(s_x[0], s_x[1], n_s_x)

        # Read the seventh record
        f.seek(6 * record_length)
        s_y = np.fromfile(f, np.float64, 2)
        s_y = np.linspace(s_y[0], s_y[1], n_s_y)

    # Read the eighth record
    f.seek(7 * record_length)
    s_z = np.fromfile(f, np.float32, n_s_z)

    # Read the ninth record
    f.seek(8 * record_length)
    r_z = np.fromfile(f, np.float32, n_r_z)

    # Read the tenth record
    f.seek(9 * record_length)
    r_r = np.fromfile(f, np.float64, n_r_r)

    # Read the pressure
    pressure = np.zeros((n_s_x, n_s_y, n_r_theta, n_s_z, n_r_z, n_r_r), dtype=np.complex64)

    for i in range(n_s_x):
        for j in range(n_s_y):
            for k in range(n_r_theta):
                for l in range(n_s_z):
                    for m in range(n_r_z):

                        # Determine the record index
                        record_index = 10 + i * n_s_y * n_r_theta * n_s_z * n_r_z + \
                                                    j * n_r_theta * n_s_z * n_r_z + \
                                                                k * n_s_z * n_r_z + \
                                                                        l * n_r_z + \
                                                                                m

                        # Read the record
                        f.seek(record_index * record_length)
                        p = np.fromfile(f, np.float32, 2 * n_r_r)

                        # Convert to a complex number
                        index = np.arange(0, 2 * n_r_r, 2)
                        pressure[i, j, k, :] = p[index] + 1j * p[index + 1]

    # Close the file
    f.close()

    # Convert distances from [km] to [m]
    s_x *= 1E3
    s_y *= 1E3
    r_r *= 1E3

    # Convert angles from [deg] to [rad]
    r_theta = np.radians(r_theta)

    return title, frequency, pressure, s_x, s_y, s_z, r_z, r_r, r_theta
