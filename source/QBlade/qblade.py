""" 
Author:   T. Moreira da Fonte Fonseca Teles
Email:    tmoreiradafont@tudelft.nl
Date:     2025-07-14
License:  GNU GPL 3.0

Manage the QBlade dynamic-link library.

Classes:
    QBlade

Functions:
    None

Exceptions:
    None
"""

from ctypes import CDLL, c_bool, c_char_p, c_double, c_int, c_void_p, POINTER
from typing import Any, Dict


class QBlade:
    """
    A class to interact with the QBlade dynamic-link library.

    Attributes:
        path : str -- path to the dynamic-link library
        lib : ctypes.CDLL -- dynamic-link library
        functions : Dict -- functions in the dynamic-link library 

    Methods:
        __init__ : initialise the QBlade class
        load_library : load the QBlade library
        unload_library : unload the QBlade library
    """

    def __init__(self, path):
        """
        Initialise the QBlade class.

        Arguments:
            path : str -- path to the dynamic-link library

        Returns:
            None
        """

        self.path = path
        self.lib = None

        # Define the function names, argument types, and return type
        self.functions: Dict[str, Dict[str, Any]] = {
            "createInstance": {"argtypes": [c_int, c_int], "restype": c_bool},
            "closeInstance": {"argtypes": None, "restype": c_void_p},
            "loadProject": {"argtypes": [c_char_p], "restype": c_void_p},
            "loadSimDefinition": {"argtypes": [c_char_p], "restype": c_void_p},
            "setOmpNumThreads": {"argtypes": [c_int], "restype": c_void_p},
            "getCustomData_at_num": {"argtypes": [c_char_p, c_double, c_int], "restype": c_double},
            "getCustomSimulationTimeData": {"argtypes": [c_char_p], "restype": c_double},
            "getWindspeed": {"argtypes": [c_double, c_double, c_double, POINTER(c_double * 3)], "restype": c_void_p},
            "getWindspeedArray": {"argtypes": [POINTER(c_double), POINTER(c_double), POINTER(c_double),POINTER(c_double), POINTER(c_double), POINTER(c_double), c_int],"restype": c_void_p,},
            "storeProject": {"argtypes": [c_char_p], "restype": c_void_p},
            "exportResults": {"argtypes": [c_int, c_char_p, c_char_p, c_char_p], "restype": c_void_p},
            "setLibraryPath": {"argtypes": [c_char_p], "restype": c_void_p},
            "setLogFile": {"argtypes": [c_char_p], "restype": c_void_p},
            "addTurbulentWind": {"argtypes": [c_double, c_double, c_double, c_double, c_int, c_double,c_double, c_char_p, c_char_p, c_int, c_double, c_double, c_bool,],"restype": c_void_p,},
            "setExternalAction": {"argtypes": [c_char_p, c_char_p, c_double, c_double, c_char_p, c_bool, c_int],"restype": c_void_p,},
            "setMooringStiffness": {"argtypes": [c_double, c_double, c_int, c_int], "restype": c_void_p},
            "loadTurbulentWindBinary": {"argtypes": [c_char_p], "restype": c_void_p},
            "setTimestepSize": {"argtypes": [c_double], "restype": c_void_p},
            "setInitialConditions_at_num": {"argtypes": [c_double, c_double, c_double, c_double, c_int],"restype": c_void_p,},
            "setRPMPrescribeType_at_num": {"argtypes": [c_int, c_int], "restype": c_void_p},
            "setRPM_at_num": {"argtypes": [c_double, c_int], "restype": c_void_p},
            "setRampupTime": {"argtypes": [c_double], "restype": c_void_p},
            "setTurbinePosition_at_num": {"argtypes": [c_double, c_double, c_double, c_double, c_double, c_double, c_int],"restype": c_void_p,},
            "getTowerBottomLoads_at_num": {"argtypes": [POINTER(c_double * 6), c_int], "restype": c_void_p},
            "initializeSimulation": {"argtypes": None, "restype": c_void_p},
            "advanceTurbineSimulation": {"argtypes": None, "restype": c_bool},
            "advanceController_at_num": {"argtypes": [POINTER(c_double * 5), c_int], "restype": c_void_p},
            "setDebugInfo": {"argtypes": [c_bool], "restype": c_void_p},
            "setUseOpenCl": {"argtypes": [c_bool], "restype": c_void_p},
            "setGranularDebug": {"argtypes": [c_bool, c_bool, c_bool, c_bool, c_bool], "restype": c_void_p},
            "setControlVars_at_num": {"argtypes": [POINTER(c_double * 5), c_int], "restype": c_void_p},
            "getTurbineOperation_at_num": {"argtypes": [POINTER(c_double * 41), c_int], "restype": c_void_p},
            "setPowerLawWind": {"argtypes": [c_double, c_double, c_double, c_double, c_double], "restype": c_void_p},
            "runFullSimulation": {"argtypes": None, "restype": c_bool},
            "setAutoCleanup": {"argtypes": [c_bool], "restype": c_void_p},
        }

        # Load the QBlade library
        self.load_library()

    def load_library(self):
        """
        Load the QBlade library.

        Arguments:
            None

        Returns:
            None
        """

        # Load the dynamic-link library
        self.lib = CDLL(self.path)

        # Bind the functions
        for name, signature in self.functions.items():

            # Get the function from the library
            function = getattr(self.lib, name)

            # Set the argument types
            function.argtypes = signature["argtypes"]

            # Set the return type
            function.restype = signature["restype"]

            # Assign the function to the class
            setattr(self, name, function)

        # Set the library path
        self.setLibraryPath(self.path.encode("utf-8"))

    def unload_library(self):
        """
        Unload the QBlade library.

        Arguments:
            None

        Returns:
            None
        """

        if self.lib:

            # Close the QBlade instance
            self.closeInstance()

            # Delete the QBlade instance
            self.lib = None
