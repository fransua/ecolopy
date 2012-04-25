#!/usr/bin/python
"""
12 Aug 2011

test functionality of abundance class
"""

__author__  = "Francois-Jose Serra"
__email__   = "francois@barrabin.org"
__licence__ = "GPLv3"
__version__ = "0.13"

import unittest
from  dataset import sample_abundance
from ecolopy_dev import Community

class Test_ewens(unittest.TestCase):
    def test_optimization(self):
        com = Community(sample_abundance)
        com.fit_model('ewens')
        ewens_model = com.get_model('ewens')
        print '* fitted to ewens model:'
        print ewens_model
        self.assertEqual(29.257363187464055, ewens_model.theta)
        self.assertEqual(107.69183285311374, ewens_model.lnL)
        self.assertEqual(14.720795937459272, ewens_model.I)
