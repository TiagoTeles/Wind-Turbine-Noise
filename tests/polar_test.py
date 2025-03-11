""" 
Author:   T. Moreira da Fonte Fonseca Teles
Email:    tmoreiradafont@tudelft.nl
Date:     2025-03-11
License:  GNU GPL 3.0

Test the Polar class.

Classes:
    TestPolar

Functions:
    None

Exceptions:
    None
"""

import os
import sys
import unittest

import numpy as np
import pandas as pd

sys.path.append(os.path.dirname(sys.path[0]))

from source.QBlade.polar import Polar


class TestPolar(unittest.TestCase):
    """
    A class to test the Polar class.

    Methods:
        test_read -- test the read method

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
        POLAR_PATH = "tests\\data\\QBlade\\NACA2412.plr"
        REFERENCE_PATH = "tests\\data\\QBlade\\polar.csv"

        # Reference attributes
        POLAR_DICT = {
            "POLARNAME": "NACA_2412_Re1.000_M0.00_N9.0 360 M",
            "FOILNAME": "NACA_2412",
            "THICKNESS": 0.12,
            "ISDECOMPOSED": True,
            "REYNOLDS": None,
        }

        # Reference polar
        reference = pd.read_csv(REFERENCE_PATH, index_col=0)

        # Create a Polar instance
        polar = Polar(POLAR_PATH)

        # Run tests
        for key, value in POLAR_DICT.items():
            self.assertEqual(polar.attributes[key], value, \
                             f"Attribute {key} does not match reference {key}.")

        pd.testing.assert_frame_equal(polar.data, reference)
