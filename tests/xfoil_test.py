""" 
Author:   T. Moreira da Fonte Fonseca Teles
Email:    tmoreiradafont@tudelft.nl
Date:     2025-03-05
License:  GNU GPL 3.0

Test the XFoil class.

Classes:
    TestXFoil

Functions:
    None

Exceptions:
    None
"""

import filecmp
import os
import sys
import unittest

import pandas as pd
import numpy as np

sys.path.append(os.path.dirname(sys.path[0]))

from source.XFOIL.xfoil import XFoil


class TestXFoil(unittest.TestCase):
    """
    A class to test the XFoil class.

    Methods:
        test_interpolate -- test the interpolate method
        test_run -- test the run method

    Attributes:
        None
    """

    def test_interpolate(self):
        """
        Test the interpolate method.

        Arguments:
            None

        Returns:
            None
        """

        # File paths
        XFOIL_PATH = "bin\\XFoil\\xfoil.exe"
        AIRFOIL_0_PATH = "tests\\data\\XFOIL\\airfoils\\NACA2412.dat"
        AIRFOIL_1_PATH = "tests\\data\\XFOIL\\airfoils\\NACA4412.dat"
        AIRFOIL_OUT_PATH = "tests\\data\\XFOIL\\airfoils\\test_airfoil.dat"
        AIRFOIL_REF_PATH = "tests\\data\\XFOIL\\reference\\interpolated_airfoil.dat"

        # Interpolate airfoils
        xfoil = XFoil(XFOIL_PATH)
        xfoil.interpolate(AIRFOIL_0_PATH, AIRFOIL_1_PATH, AIRFOIL_OUT_PATH, 0.5)

        # Run test
        self.assertTrue(filecmp.cmp(AIRFOIL_OUT_PATH, AIRFOIL_REF_PATH, shallow=False), \
                        "Interpolated airfoil does not match reference airfoil!")

        # Remove output file
        os.remove(AIRFOIL_OUT_PATH)

    def test_run(self):
        """
        Test the run method.

        Arguments:
            None
        
        Returns:
            None
        """

        # File paths
        XFOIL_PATH = "bin\\XFoil\\xfoil.exe"
        AIRFOIL_0_PATH = "tests\\data\\XFOIL\\airfoils\\NACA2412.dat"
        BL_TOP_REF_PATH = "tests\\data\\XFOIL\\reference\\boundary_layer_top.csv"
        BL_BOT_REF_PATH = "tests\\data\\XFOIL\\reference\\boundary_layer_bot.csv"

        # Load reference output
        reference_top = pd.read_csv(BL_TOP_REF_PATH, index_col=0)
        reference_bot = pd.read_csv(BL_BOT_REF_PATH, index_col=0)

        # Run XFOIL
        xfoil = XFoil(XFOIL_PATH)
        actual_top, actual_bot = xfoil.run(AIRFOIL_0_PATH, 1E6, 0.2, np.radians(10), it=100)

        # Run test
        self.assertTrue(actual_top.equals(reference_top), \
                        "XFOIL output does not match reference output!")

        self.assertTrue(actual_bot.equals(reference_bot), \
                        "XFOIL output does not match reference output!")
