import sys

def parse(file, key, index, arg_type):
    """ 
    Searches for a key in a file and returns the value at a specified index.

    Arguments:
        file -- the file to search in
        key -- the string to search for
        index -- the index of the value to return
        arg_type -- the type of the value to return

    Returns:
        the value at the specified index with a specified key and type
    """

    for line in file:
        if key in line:
            if arg_type == str:
                return str(line.split()[index])
            elif arg_type == int:
                return int(line.split()[index])
            elif arg_type == float:
                return float(line.split()[index])
            elif arg_type == bool:
                return line.split()[index] == "true"
            else:
                print("Invalid argument type!")
                sys.exit(1)

    print("Key not found!")
    sys.exit(1)
