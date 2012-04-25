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
from dataset     import sample_abundance
from ecolopy_dev import Community

class Test_etienne(unittest.TestCase):
    def test_get_kda(self):
        pass
        
    def test_optimization(self):
        com = Community(sample_abundance)
        com.fit_model('etienne')
        etienne_model = com.get_model('etienne')
        print '* fitted to etienne model:'
        print etienne_model
        self.assertEqual(round(66.1478972798     , 5), round(etienne_model.theta, 5))
        self.assertEqual(round(98.171564034001662, 5), round(etienne_model.lnL  , 5))
        self.assertEqual(round(134.10853894914644, 2), round(etienne_model.I    , 2))

    def test_randomization(self):
        com = Community(sample_abundance)       
        com.fit_model('etienne')                
        etienne_model = com.get_model('etienne')
        self.assertEqual(sum(etienne_model.random_community()), com.J)
        self.assertEqual(sum(etienne_model.random_community(2500)), 2500)