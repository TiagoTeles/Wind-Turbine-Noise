""" 
Author:   T. Moreira da Fonte Fonseca Teles
Email:    tmoreiradafont@tudelft.nl
Date:     2025-03-05
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

    def test_E(self):
        """
        Test the E function.

        Arguments:
            None
        
        Returns:
            None
        """
