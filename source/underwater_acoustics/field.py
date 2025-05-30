""" 
Author:   T. Moreira da Fonte Fonseca Teles
Email:    tmoreiradafont@tudelft.nl
Date:     2025-05-30
License:  GNU GPL 3.0

Write .flp files.

Classes:
    None

Functions:
    write_flp
    uniform_mesh

Exceptions:
    None
"""

import numpy as np


def write_flp(path, title, max_modes, s_x, s_y, s_z, r_z, r_r, r_theta, nodes, elements):
    """
    Write the field file for a Nx2D simulation.

    Arguments:
        path : str -- path to the field file
        title : str -- title of the field file
        max_modes : int -- maximum number of modes to be computed, [-]
        s_x : np.ndarray -- source x coordinates, [m]
        s_y : np.ndarray -- source y coordinates, [m]
        s_z : np.ndarray -- source z coordinates, [m]
        r_z : np.ndarray -- receiver z coordinates, [m]
        r_r : np.ndarray -- receiver ranges, [m]
        r_theta : np.ndarray -- receiver bearings, [rad]
        nodes : pd.DataFrame -- nodes
        elements : np.ndarray -- elements
    """

    # Convert positions from [km] to [m]
    s_x /= 1E3
    s_y /= 1E3
    r_r /= 1E3
    nodes["x"] /= 1E3
    nodes["y"] /= 1E3

    # Convert angles from [rad] to [deg]
    r_theta = np.degrees(r_theta)

    # Create the field file
    f = open(path, "w", encoding="utf-8")

    # Write the title
    f.write(f"{title}\n")

    # Write the Field3D routine
    f.write("STD\n")

    # Write the maximum number of modes
    f.write(f"{max_modes}\n")

    # Write the source x coordinates
    f.write(f"{len(s_x)}\n")
    f.write(f"{s_x[0]} {s_x[-1]} /\n")

    # Write the source y coordinates
    f.write(f"{len(s_y)}\n")
    f.write(f"{s_y[0]} {s_y[-1]} /\n")

    # Write the source z coordinates
    f.write(f"{len(s_z)}\n")
    f.write(f"{s_z[0]} {s_z[-1]} /\n")

    # Write the receiver z coordinates
    f.write(f"{len(r_z)}\n")
    f.write(f"{r_z[0]} {r_z[-1]} /\n")

    # Write the receiver ranges
    f.write(f"{len(r_r)}\n")
    f.write(f"{r_r[0]} {r_r[-1]} /\n")

    # Write the receiver bearings
    f.write(f"{len(r_theta)}\n")
    f.write(f"{r_theta[0]} {r_theta[-1]} /\n")

    # Write the nodes
    n_nodes = nodes.shape[0]
    f.write(f"{n_nodes}\n")

    for i in range(n_nodes):

        node_x = nodes["x"].iloc[i]
        node_y = nodes["y"].iloc[i]
        node_filename = nodes["mod_file"].iloc[i]

        f.write(f"{node_x} {node_y} '{node_filename}'\n")

    # Write the elements
    n_elements = elements.shape[0]
    f.write(f"{n_elements}\n")

    for i in range(n_elements):

        node_1 = elements[i, 0] + 1
        node_2 = elements[i, 1] + 1
        node_3 = elements[i, 2] + 1

        f.write(f"{node_1} {node_2} {node_3}\n")

    # Close the file
    f.close()


def uniform_mesh(N_x, N_y):
    """
    Determine the elements of the uniform triangular mesh.

    Arguments:
        N_x : int -- number of nodes in the x-direction, [-]
        N_y : int -- number of nodes in the y-direction, [-]

    Returns:
        elements : np.ndarray -- mesh elements
    """

    elements = np.zeros((2 * (N_x - 1) * (N_y - 1), 3), dtype=int)

    for i in range(N_x-1):
        for j in range(N_y-1):

            # Node indices
            sw_node = j * N_x + i
            se_node = j * N_x + (i + 1)
            nw_node = (j + 1) * N_x + i
            ne_node = (j + 1) * N_x + (i + 1)

            # Element indices
            lower_element = 2 * (j * (N_x - 1) + i)
            upper_element = 2 * (j * (N_x - 1) + i) + 1

            # Lower element
            elements[lower_element, 0] = sw_node
            elements[lower_element, 1] = se_node
            elements[lower_element, 2] = ne_node

            # Upper element
            elements[upper_element, 0] = sw_node
            elements[upper_element, 1] = ne_node
            elements[upper_element, 2] = nw_node

    return elements
