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
        
    def test_optimization(self):
        com = Community(sample_abundance)
        com.fit_model('lognormal')
        model = com.get_model('lognormal')
        print '* fitted to lognormal model:'
        print model
        self.assertEqual(round(2.14539 , 5), round(model.mu, 5))
        self.assertEqual(round(1.28561 , 5), round(model.sd  , 5))

    def test_randomization(self):
        com = Community(sample_abundance)       
        com.fit_model('lognormal')                
        model = com.get_model('lognormal')
        self.assertTrue(sum(model.random_community())>1000)
    