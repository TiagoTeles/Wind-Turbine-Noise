""" 
Author:   T. Moreira da Fonte Fonseca Teles
Email:    tmoreiradafont@tudelft.nl
Date:     2025-05-28
License:  GNU GPL 3.0

Setup the .flp file.

Classes:
    None

Functions:
    read_shade

Exceptions:
    None
"""

import numpy as np

def read_shade(filename):
    """
    Read .shd files.

    Arguments:
        filename : str -- path to the .shd file

    Returns:
        pressure : np.ndarray -- pressure, [Pa]
        frequency : np.ndarray -- frequency, [Hz]
        theta : np.ndarray -- bearing angle, [rad]
        source_position : tuple -- source positions, [m]
        receiver_position : tuple -- receiver positions, [m]
    """

    # Open the file
    f = open(filename, "rb")

    # Read the record length in bytes
    record_length = 4 * np.fromfile(f, np.int32, 1)[0]

    # Read the plot title
    title = f.read(80)

    # Read the first record
    f.seek(1 * record_length)
    plot_type = f.read(10)

    # Read the second record
    f.seek(2 * record_length)
    n_frequencies = np.fromfile(f, np.int32, 1)[0]
    n_theta = np.fromfile(f, np.int32, 1)[0]
    n_s_x = np.fromfile(f, np.int32, 1)[0]
    n_s_y = np.fromfile(f, np.int32, 1)[0]
    n_s_z = np.fromfile(f, np.int32, 1)[0]
    n_r_z = np.fromfile(f, np.int32, 1)[0]
    n_r_r = np.fromfile(f, np.int32, 1)[0]
    frequency_0 = np.fromfile(f, np.float64, 1)[0]
    attenuation = np.fromfile(f, np.float64, 1)[0]

    # Read the third record
    f.seek(3 * record_length)
    frequency = np.fromfile(f, np.float64, n_frequencies)

    # Read the fourth record
    f.seek(4 * record_length)
    theta = np.fromfile(f, np.float64, n_theta)

    if plot_type[0 : 1] != "TL":

        # Read the fifth record
        f.seek(5 * record_length)
        x_s = np.fromfile(f, np.float64, n_s_x)

        # Read the sixth record
        f.seek(6 * record_length)
        y_s = np.fromfile(f, np.float64, n_s_y)

    else:

        # Read the fifth record
        f.seek(5 * record_length)
        x_s = np.fromfile(f, np.float64, 2)
        x_s = np.linspace(x_s[0], x_s[1], n_s_x)

        # Read the sixth record
        f.seek(6 * record_length)
        y_s = np.fromfile(f, np.float64, 2)
        y_s = np.linspace(y_s[0], y_s[1], n_s_y)

    # Read the seventh record
    f.seek(7 * record_length)
    z_s = np.fromfile(f, np.float32, n_s_z)

    # Read the eighth record
    f.seek(8 * record_length)
    z_r = np.fromfile(f, np.float32, n_r_z)

    # Read the ninth record
    f.seek(9 * record_length)
    r_r = np.fromfile(f, np.float64, n_r_r)

    # Read the pressure
    pressure = np.zeros((n_theta, n_s_z, n_r_z, n_r_r), dtype=np.complex64)

    # Read the pressure
    for i in range(n_theta):
        for j in range(n_s_z):
            for k in range(n_r_z):

                record_index = 10  + i * n_s_z * n_r_z \
                                           + j * n_r_z \
                                                   + k

                # Read the record
                f.seek(record_index * record_length)

                p = np.fromfile(f, np.float32, 2 * n_r_r)
                index = np.arange(0, 2 * n_r_r, 2)

                pressure[i, j, k, :] = p[index] + 1j * p[index + 1]

    # Close the file
    f.close()

    # TODO: CONVERT UNITS OF ANGLES AND DISTANCES

    return pressure, frequency, theta, (x_s, y_s, z_s), (z_r, r_r)
