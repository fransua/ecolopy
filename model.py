#!/usr/bin/python
"""
13 Aug 2011


"""

__author__  = "Francois-Jose Serra"
__email__   = "francois@barrabin.org"
__licence__ = "GPLv3"
__version__ = "0.0"


class Model (object):
    def __init__(self, name):
        self.name = name

    def get_likelihood(self):
        pass

    def get_optimal_params(self):
        pass

    
