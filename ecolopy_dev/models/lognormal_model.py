#!/usr/bin/python
"""
21 Sep 2011


"""

__author__  = "Francois-Jose Serra"
__email__   = "francois@barrabin.org"
__licence__ = "GPLv3"
__version__ = "0.13"


from ecolopy_dev.models import EcologicalModel
from ecolopy_dev.utils  import table, lpoch
from math               import log, lgamma
from numpy              import mean, std
from random             import lognormvariate

class LognormalModel(EcologicalModel):
    def __init__(self, community, **kwargs):
        super(LognormalModel, self).__init__(community, **kwargs)
        self._parameters['sd'] = None
        self._parameters['mu'] = None        
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
        self._parameters['mu'] = mean([float (log (x)) for x in self.community.abund])
        self._parameters['sd'] = std ([float (log (x)) for x in self.community.abund])
        self._lnL              = None


    def random_community(self, inds=None, mu=None, sd=None):
        '''
        generates random lognormal distribution

        :argument inds: number of individuals in community (J)
        :argument sd: usually standard deviation of the distribution of abundances
        :argument mu: usually mean of the distribution of abundances

        :returns: distribution of abundance (list)
        '''
        inds = inds or self.community.S
        mu = mu or self._parameters['mu']
        sd = sd or self._parameters['sd']
        return [lognormvariate (mu, sd) for _ in xrange(inds)]
        
    def get_mu(self):
        return self._parameters['mu']

    mu = property(get_mu, doc="mu")

    def get_sd(self):
        return self._parameters['sd']

    sd = property(get_sd, doc="standard deviation")

        