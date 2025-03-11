""" 
Author:   T. Moreira da Fonte Fonseca Teles
Email:    tmoreiradafont@tudelft.nl
Date:     2025-03-11
License:  GNU GPL 3.0

Test miscellaneous functions.

Classes:
    TestMisc

Functions:
    None

Exceptions:
    None
"""

import os
import sys
import unittest

import numpy as np

sys.path.append(os.path.dirname(sys.path[0]))

from source.misc import octave, E


class TestMisc(unittest.TestCase):
    """
    A class to test miscellaneous functions.

    Methods:
        test_octave -- test the octave function
        test_E -- test the E function

    Attributes:
        None
    """

    def test_octave(self):
        """
        Test the octave function.

        Arguments:
            None

        Returns:
            None
        """

        # Reference values
        reference = np.array([19.953, 25.119, 31.623, 39.811, 50.119, 63.096, 79.433, 100.00,
                              125.89, 158.49, 199.53, 251.19, 316.23, 398.11, 501.19, 630.96,
                              794.43, 1000.0, 1258.9, 1584.9, 1995.3, 2511.9, 3162.3, 3981.1,
                              5011.9, 6309.6, 7943.30, 10000,  12589,  15849, 19953])

        # Actual values
        actual, _, _ = octave(20, 20000, 1000, True)

        np.testing.assert_allclose(actual, reference, rtol=1e-3, \
                                   err_msg="Frequencies do not match reference frequencies!")

    def test_E(self):
        """
        Test the E function.

        Arguments:
            None
        
        Returns:
            None
        """

        # X-axis values
        x = np.arange(1E-6, 10, 1E-6)

        # Reference values
        reference = np.cumsum(np.exp(-1j * x) / np.sqrt(2 * np.pi * x) * 1E-6)

        # Actual values
        actual = E(x)

        np.testing.assert_allclose(actual, reference, atol=1E-3, \
                                   err_msg="Fresnel integrals do not match reference values!")
