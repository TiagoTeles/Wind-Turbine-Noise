""" 
Author:   T. Moreira da Fonte Fonseca Teles
Email:    tmoreiradafont@tudelft.nl
Date:     2025-11-11
License:  GNU GPL 3.0

Manage the QBlade library.

Classes:
    QBlade

Functions:
    None

Exceptions:
    None
"""

from ctypes import CDLL, c_bool, c_char_p, c_int, c_void_p


class QBlade:
    """
    A class to manage the QBlade library.

    Methods:
        __init__ : initialise the QBlade class
        set_library_path : set the QBlade library path
        create_instance : create the QBlade instance
        load_sim_definition : load the QBlade simulation definition file
        initialise_simulation : initialise the QBlade simulation
        run_full_simulation : run the QBlade simulation
        export_results : export the QBlade simulation results
        close_instance : close the QBlade instance

    Attributes:
        path : str -- path to the QBlade library
        lib : ctypes.CDLL -- QBlade library
    """

    def __init__(self, path):
        """
        Initialise the QBlade class.

        Parameters:
            path : str -- path to the QBlade library

        Returns:
            None
        """

        self.path = path

        # Load the QBlade library
        self.lib = CDLL(self.path)

        # Set the QBlade library path
        self.set_library_path(self.path)

    def set_library_path(self, path):
        """
        Set the QBlade library path.

        Parameters:
            path : str -- path to the QBlade library

        Returns:
            None
        """

        # Set the argument and return types
        self.lib.setLibraryPath.argtypes = [c_char_p]
        self.lib.setLibraryPath.restype = c_void_p

        # Call the function
        self.lib.setLibraryPath(path.encode("utf-8"))

        return None

    def create_instance(self, cl_device, group_size):
        """
        Create the QBlade instance.

        Parameters:
            cl_device : int -- OpenCL device
            group_size : int -- OpenCL work-group size

        Returns:
            success : bool -- status
        """

        # Set the argument and return types
        self.lib.createInstance.argtypes = [c_int, c_int]
        self.lib.createInstance.restype = c_bool

        # Call the function
        success = self.lib.createInstance(cl_device, group_size)

        return success

    def load_sim_definition(self, path):
        """
        Load the QBlade simulation definition file.

        Parameters:
            path : str -- path to the QBlade simulation definition file

        Returns:
            None
        """

        # Set the argument and return types
        self.lib.loadSimDefinition.argtypes = [c_char_p]
        self.lib.loadSimDefinition.restype = c_void_p

        # Call the function
        self.lib.loadSimDefinition(path.encode("utf-8"))

        return None

    def initialise_simulation(self):
        """
        Initialise the QBlade simulation.

        Parameters:
            None

        Returns:
            None
        """

        # Set the argument and return types
        self.lib.initializeSimulation.argtypes = None
        self.lib.initializeSimulation.restype = c_void_p

        # Call the function
        self.lib.initializeSimulation()

        return None

    def run_full_simulation(self):
        """
        Run the QBlade simulation.

        Parameters:
            None

        Returns:
            success : bool -- status
        """

        # Set the argument and return types
        self.lib.runFullSimulation.argtypes = None
        self.lib.runFullSimulation.restype = c_bool

        # Call the function
        success = self.lib.runFullSimulation()

        return success

    def export_results(self, type, directory, name, filter):
        """
        Export the QBlade simulation results.

        Parameters:
            type : int -- file type
            directory : str -- directory
            name : str -- name
            filter : str -- filter

        Returns:
            None
        """

        # Set the argument and return types
        self.lib.exportResults.argtypes = [c_int, c_char_p, c_char_p, c_char_p]
        self.lib.exportResults.restype = c_void_p

        # Call the function
        self.lib.exportResults(type, directory.encode("utf-8"), name.encode("utf-8"), \
                               filter.encode("utf-8"))

        return None

    def close_instance(self):
        """
        Close the QBlade instance.

        Parameters:
            None

        Returns:
            None
        """

        # Set the argument and return types
        self.lib.closeInstance.argtypes = None
        self.lib.closeInstance.restype = c_void_p

        # Call the function
        self.lib.closeInstance()

        # Delete the QBlade instance
        self.lib = None

        return None
