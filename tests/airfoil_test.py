""" 
Author:   T. Moreira da Fonte Fonseca Teles
Email:    tmoreiradafont@tudelft.nl
Date:     2025-03-11
License:  GNU GPL 3.0

Test the Airfoil class.

Classes:
    TestAirfoil

Functions:
    None

Exceptions:
    None
"""

import os
import sys
import unittest

import pandas as pd

sys.path.append(os.path.dirname(sys.path[0]))

from source.QBlade.airfoil import Airfoil


class TestAirfoil(unittest.TestCase):
    """
    A class to test the Airfoil class.

    Methods:
        test_read -- test the read method
        test_thickness -- test the thickness method

    Attributes:
        None
    """

    def test_read(self):
        """
        Test the read method.

        Arguments:
            None

        Returns:
            None
        """

        # File path
        AIRFOIL_PATH = "tests\\data\\QBlade\\NACA2412.afl"
        REFERENCE_PATH = "tests\\data\\QBlade\\coordinates.csv"

        # Reference coordinates
        reference = pd.read_csv(REFERENCE_PATH, index_col=0)

        # Create an Airfoil instance
        airfoil = Airfoil(AIRFOIL_PATH)

        # Run tests
        self.assertEqual(airfoil.attributes["AIRFOILNAME"], "NACA 2412", \
                         "Airfoil name does not match reference name.")

        pd.testing.assert_frame_equal(airfoil.data, reference)

    def test_thickness(self):
        """"
        Test the thickness method.

        Arguments:
            None

        Returns:
            None
        """

        # File path
        AIRFOIL_PATH = "tests\\data\\QBlade\\NACA2412.afl"

        # Reference thickness
        tc_01 = 0.034017165
        tc_10 = 0.093615533

        # Create an Airfoil instance
        airfoil = Airfoil(AIRFOIL_PATH)

        # Run tests
        self.assertAlmostEqual(airfoil.thickness(0.01), tc_01, \
                               msg="Thickness does not match reference thickness!")

        self.assertAlmostEqual(airfoil.thickness(0.10), tc_10, \
                               msg="Thickness does not match reference thickness!")
