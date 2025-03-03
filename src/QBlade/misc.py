"""
Author:   T. Moreira da Fonte Fonseca Teles
Email:    tmoreiradafont@tudelft.nl
Date:     2025-02-26
License:  GNU GPL 3.0

Helper functions for reading and writing files.

Classes:
    None

Functions:
    read
    write

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
        arg_type : dict -- argument type

    Returns:
        value : Any -- value of the key
    """

    # Go to file start
    file.seek(0)

    # Search for the key
    for line in file:
        if re.search(rf'\b{re.escape(key)}\b', line) is not None:

            value = line[0:line.find(key)].strip()

            if arg_type in {str, int, float}:
                return arg_type(value)

            elif arg_type == bool:
                return value == "true"

            elif arg_type is None:
                return None

            else:
                print("Invalid argument type!")
                sys.exit(1)

    print("Key not found!")
    sys.exit(1)

def write(file, key, value):
    """
    Set the value of a key in a file.

    Arguments:
        file : file -- file object
        key : str -- key to search for
        value : Any -- value to set

    Returns:
        None
    """

    # Go to file start
    file.seek(0)

    # Search for the key
    lines = []

    new_value = None

    for line in file:
        if re.search(rf'\b{re.escape(key)}\b', line) is not None:

            # Determine the previous value
            old_value = line[0:line.find(key)].strip()

            # Determine the new value
            if type(value) in {str, int, float}:
                new_value = str(value)

            elif isinstance(value, bool):
                if value:
                    new_value = "true"
                else:
                    new_value = "false"

            elif value is None:
                new_value = old_value

            else:
                print("Invalid argument type!")
                sys.exit(1)

            # Add padding
            if len(old_value) > len(new_value):
                new_value = new_value + " " * (len(old_value) - len(new_value))
            else:
                old_value = old_value + " " * (len(new_value) - len(old_value))

            # Replace the old value with the new value
            lines.append(line.replace(old_value, new_value, 1))

        else:
            lines.append(line)

    # Write the new values to the file
    file.seek(0)
    file.truncate(0)
    file.writelines(lines)

    if new_value is None:
        print("Key not found!")
        sys.exit(1)
