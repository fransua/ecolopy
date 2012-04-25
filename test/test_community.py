#!/usr/bin/python
"""
04 Jul 2011


"""

__author__  = "Francois-Jose Serra"
__email__   = "francois@barrabin.org"
__licence__ = "GPLv3"
__version__ = "0.13"

import unittest
from dataset import sample_abundance

from ecolopy_dev import Community

class TestCommunity(unittest.TestCase):
    """
    Test Community basics
    """
    def test_reading_file(self):
        com = Community(sample_abundance)
        self.assertEqual(com.S, 127)
        self.assertEqual(com.J, 2179)
        print '* abundance loaded into community:'
        print com


if __name__ == '__main__':
    unittest.main()
