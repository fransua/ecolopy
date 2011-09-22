#!/usr/bin/python
"""
21 Sep 2011


"""

__author__  = "Francois-Jose Serra"
__email__   = "francois@barrabin.org"
__licence__ = "GPLv3"
__version__ = "0.13"


from random_neutral import rand_neutral_ewens
from random_neutral import rand_neutral_etienne
from random_neutral import rand_lognormal

class EcologicalModel(object):
    """
    """
    def __init__(self, name, **kwargs):
        self.name  = name
        for key in kwargs:
            setattr (self, key, kwargs[key])
        if not 'theta' in kwargs or not 'I' in kwargs:
            raise Exception ('Must supply theta and I values\n')


    def rand_neutral (self, inds):
        """
        """
        if self.name == 'ewens':
            return rand_neutral_ewens (inds, self.theta)
        elif self.name == 'etienne':
            return rand_neutral_etienne (inds, self.theta, self.I)
        elif self.name == 'lognorm':
            return rand_lognormal (inds, self.theta, self.I)
        else:
            raise Exception ('No random simulation available for this model\n')
        

