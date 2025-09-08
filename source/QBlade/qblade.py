""" 
Author:   T. Moreira da Fonte Fonseca Teles
Email:    tmoreiradafont@tudelft.nl
Date:     2025-08-22
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
    A class to interact with the QBlade library.

    Methods:
        __init__ : initialise the QBlade class
        set_library_path : set the QBlade library path
        create_instance : create the QBlade instance
        load_sim_definition : load the simulation definition file
        initialise_simulation : initialise the simulation
        run_full_simulation : run the simulation
        export_results : export the simulation results
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

        # Load the library
        self.lib = CDLL(self.path)

        # Set the library path
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
            group_size : int -- work group size

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
        Load the simulation definition file.

        Parameters:
            path : str -- path to the simulation definition file

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
        Initialise the simulation.

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
        Run the simulation.

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
        Export the simulation results.

        Parameters:
            type : int -- results file type
            directory : str -- results directory
            name : str -- results file name
            filter : str -- results filter

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
