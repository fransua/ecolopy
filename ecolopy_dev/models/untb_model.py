#!/usr/bin/python
"""
21 Sep 2011


"""

__author__  = "Francois-Jose Serra"
__email__   = "francois@barrabin.org"
__licence__ = "GPLv3"
__version__ = "0.13"


from ecolopy_dev.models import EcologicalModel
from ecolopy_dev.utils  import table
from math               import lgamma, log

class UNTBModel(EcologicalModel):
    """
    Class representing Ecological models

    :argument name: name of the class, can be either ewens, etienne or lognorm
    :returns: EcologicalModel object
    """
    def __init__(self, community, **kwargs):
        super(UNTBModel, self).__init__(community, **kwargs)
        self.__compute_factor()
        if not 'kda' in kwargs:
            self._kda = None
        self._parameters['theta']  = None
        self._parameters['I']      = None
        self._parameters['m']      = None
        self._factor = None
        self.__compute_factor()
        
    def __compute_factor(self):
        self._factor = lgamma (self.community.J + 1)
        phi = table(self.community.abund)
        phi += [0] * int (max (self.community.abund) - len (phi))
        for spe in xrange (self.community.S):
            self._factor -= log (max (1, self.community.abund[spe]))
        for spe in xrange (int(max(self.community.abund))):
            self._factor -= lgamma (phi[spe] + 1)

            
    def get_theta(self):
        return self._parameters['theta']

    theta = property(get_theta, doc="Fundamental biodiversity number")

    def get_I(self):
        return self._parameters['I']

    I = property(get_I, doc="Immigration")

    def get_m(self):
        return self._parameters['m']
        
    m = property(get_m, doc="Mutation rate")
