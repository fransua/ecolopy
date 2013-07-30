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
        self.assertEqual(round(29.257363187464055, 4), round(ewens_model.theta, 4))
        self.assertEqual(round(107.69183285311374, 4), round(ewens_model.lnL, 4))
        self.assertEqual(14.720795937459272, ewens_model.I)
