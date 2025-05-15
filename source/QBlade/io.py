"""
Author:   T. Moreira da Fonte Fonseca Teles
Email:    tmoreiradafont@tudelft.nl
Date:     2025-05-15
License:  GNU GPL 3.0

Functions for reading QBlade files.

Classes:
    None

Functions:
    read

Exceptions:
    None
"""

import re
import sys


def read(file, key, arg_type):
    """
    Get the value of a key in a file.

    Arguments:
        file : file -- file object
        key : str -- key to search for
        arg_type : type -- argument type

    Returns:
        value : Any -- value of the key
    """

    # Go to the beginning of the file
    file.seek(0)

    # Search for the key
    for line in file:
        if re.search(rf'\b{re.escape(key)}\b', line) is not None:

            value = line[0:line.find(key)].strip()

            if arg_type in {str, int, float}:
                return arg_type(value)

            elif arg_type in {bool}:
                return value == "true"

            else:
                print("Invalid argument type!")
                sys.exit(1)

    print("Key not found!")
    sys.exit(1)
