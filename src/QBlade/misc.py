"""
Author:   T. Moreira da Fonte Fonseca Teles
Email:    tmoreiradafont@tudelft.nl
Date:     2025-02-19
License:  GNU GPL 3.0

Helper functions for reading and writing files.

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
        arg_type : dict -- dictionary of keys

    Returns:
        value : Any -- value of the key
    """

    # Go to file start
    file.seek(0)

    # Search for key
    for line in file:
        if re.search(rf'\b{re.escape(key)}\b', line) is not None:

            value = line[0:line.find(key)].strip()

            if arg_type == str:
                return str(value)

            elif arg_type == int:
                return int(value)

            elif arg_type == float:
                return float(value)

            elif arg_type == bool:
                return value == "true"

            elif arg_type == None:
                print(f"Skipping key: '{key}'!")
                return None

            else:
                print("Invalid argument type!")
                sys.exit(1)

    print("Key not found!")
    sys.exit(1)


def write(file, key, arg_type, index):
    pass
