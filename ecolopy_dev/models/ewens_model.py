#!/usr/bin/python
"""
21 Sep 2011


"""

__author__  = "Francois-Jose Serra"
__email__   = "francois@barrabin.org"
__licence__ = "GPLv3"
__version__ = "0.13"


from ecolopy_dev.models import UNTBModel
from ecolopy_dev.utils  import table, lpoch
from scipy.optimize     import golden
from math               import log, lgamma
from random             import random


class EwensModel(UNTBModel):
    def __init__(self, community, **kwargs):
        super(EwensModel, self).__init__(community, **kwargs)
        self.optimize()
    
    def likelihood (self, theta):
        '''
        get likelihood value of Ewens according to Parthy/Tetame (Jabot 2008)
        
        :argument theta: value of theta
        :returns: likelihood
        '''
        factor = lgamma (self.community.J + 1)
        phi = table (self.community.abund)
        phi += [0] * int (max (self.community.abund) - len (phi))
        for spe in xrange (self.community.S):
            factor -= log (max (1, self.community.abund[spe]))
        for spe in xrange (max (self.community.abund)):
            factor -= lgamma (phi[spe] + 1)
        lnl = lpoch (theta, self.community.J) - log (theta) * self.community.S - factor
        self._factor = factor
        return lnl

    def optimize (self):
        '''
        Main function to optimize theta using theta likelihood function, according to Ewens
        model.
        '''
        theta_like   = lambda x: -self._ewens_theta_likelihood (x)
        self._parameters['theta'] = golden (theta_like, brack=[.01/self.community.J, self.community.J])
        self._parameters['m']     = self.theta / self.community.J / 2
        self._parameters['I']     = self.m * (self.community.J - 1) / (1 - self.m)
        self._lnL                = self.likelihood (self.theta)

    def _ewens_theta_likelihood (self, theta):
        '''
        returns the likelihood of theta for a given dataset
        
        '''
        if theta < 0:
            return float ('-inf')
        return self.community.S * log(theta) + lgamma(theta) - lgamma(theta + self.community.J)


    def random_community(self, inds=None, theta=None):
        '''
        generates random distribution according to J and theta
        
        :argument inds: number of individuals in community (J)
        :argument theta: corresponding to the model
    
        :returns: distribution of abundance (list)
        
        '''
        theta = float (theta) if theta else self.theta
        inds = inds if inds else self.community.J
        out = [0] * int (inds)
        out [0] = spp = 1
        for ind in xrange (inds):
            if random () < theta/(theta + ind):
                spp += 1
                out[ind] = spp
            else:
                out[ind] = out [int (random () * ind)]
        return table (out, spp)
        
        