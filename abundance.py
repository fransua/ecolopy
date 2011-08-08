#!/usr/bin/python
"""
04 Jul 2011


"""

__author__  = "Francois-Jose Serra"
__email__   = "francois@barrabin.org"
__licence__ = "GPLv3"
__version__ = "0.0"

try:
    from math import lgamma
except ImportError:
    from scipy.special import gammaln as lgamma
from numpy         import log, float128, exp
from random        import random
from scipy.stats   import chisqprob
from scipy         import optimize
from scipy.special.orthogonal import poch

from utils         import table, get_kda

class Abundance (object):
    '''
    Compute theta m I ect... for a given data sample
    '''
    def __init__ (self, data, J_tot=None):
        if type (data) != list:
            self.data_path = data
            data = self.parse_infile ()
        self.abund     = data [:]
        self.J         = sum (data)
        self.S         = len (data)
        self.theta     = self._ewens_optimal_theta ()
        self.m         = self.theta / self.J / 2
        self.I         = self.m * float (self.J - 1) / (1 - self.m)
        self.J_tot     = J_tot if J_tot else self.J * 3
        self.H         = self._shannon_entropy ()
        self.factor    = None
        self.ewens_lnl = None 
        self.K         = None

    def _lrt_ewens_etienne(self):
        '''
        NOT WORKING
        ***********
        likelihood-ratio-test between two neutral models
        returns p-value of rejection of alternative model
        '''
        chisqprob(2 * (self.ewens_likelihood() - 0), df=1)

    def ewens_likelihood (self):
        '''
        get likelihood value of Ewens according to parthy/tetame
        returns likelihood
        '''
        factor = lgamma (self.J+1)
        phi = table (self.abund)
        phi += [0] * (max (self.abund)-len (phi))
        for spe in xrange (self.S):
            factor -= log (max (1, self.abund[spe]))
        for spe in xrange (max (self.abund)):
            factor -= lgamma (phi[spe] + 1)
        self.ewens_lnl = poch (self.theta,
                               self.J) - log (self.theta) * \
                               self.S - factor
        self.factor    = factor
        return self.ewens_lnl

    def _ewens_optimal_theta (self):
        '''
        optimize theta using theta likelihood function
        '''
        theta_like = lambda x: -self._ewens_theta_likelihood (x)
        return optimize.golden(theta_like, brack=[.01/self.J, self.J])

    def _ewens_theta_likelihood (self, theta):
        '''
        returns the likelihood of theta for a given dataset
        '''
        if theta < 0:
            return float ('-inf')
        theta = float (theta)
        return self.S * log (theta) + lgamma (theta) - lgamma (theta + self.J)

    def _shannon_entropy (self):
        '''
        computes shannon entropy
        '''
        H = 0
        for spe in self.abund:
            H -= spe * log (spe)
        return (H + self.J * log(self.J)) / self.J

    def test_neutrality (self, gens=100):
        '''
        test for neutrality comparing shanon entropy
        returns p_value
        '''
        def local_shannon (data, J):
            '''
            same as self.shannon entropy
            '''
            H = 0
            for spe in data:
                H -= spe * log (spe)
            return (H + J * log(J)) / J
        pval = 0
        for _ in xrange (gens):
            pval += local_shannon (self.rand_neutral(), self.J) < self.H
        return float (pval)/gens
    
    def rand_neutral (self):
        '''
        J should be the number of individuals in the dataset, and theta,
        the optimal theta of the dataset accod
        J = 376
        theta = 9.99
        '''
        theta = float (self.theta)
        out = [0] * self.J
        out [0] = spp = 1
        for ind in xrange (1, self.J):
            if random () < theta/(theta + ind):
                spp += 1
                out [ind] = spp
            else:
                out [ind] = out [int (random () * ind)]
        return table (out, spp)
        
    def parse_infile (self):
        '''
        parse infile and return list of abundances
        '''
        abundances = []
        for line in open (self.data_path):
            abundances.append (int (line.strip()))
        return abundances

    def etienne_likelihood(self, x):
        '''
        logikelihood function where
          x[0] = theta
          x[1] = I
          K[A] = log (K(D,A))
          K[A] = K(D,A) in Etienne paper
        cf Etienne, 2005
        x = [48.1838, 1088.19]
        '''
        divisor    = 0.0
        newdivisor = 0.0
        summand0   = 0.0
        summand1   = 0.0
        lnum = float128 (11300)
        if not self.factor:
            self.ewens_likelihood()
        if not self.K:
            self.K = get_kda (self.abund)
        for A in xrange (self.J-self.S):
            lsummand = float128 (self.factor + log(x[0]) * self.S - \
                                 poch (x[1], int (self.J)) + self.K [A] + \
                                 (A + self.S) * log (x[1]) - \
                                 poch (x[0], A + self.S) - divisor)
            if lsummand > 11300:
                newdivisor = float128 (lsummand - 11300)
                divisor += float128 (newdivisor)
                summand0 = float128 (summand0 / exp (newdivisor) + exp (lnum))
            else:
                if lsummand > -11333.2 and summand0 < 11330:
                    summand0 += exp (lsummand)
                else:
                    if summand0 > 11330:
                        divisor  = float128 (1 + divisor)
                        summand0 = float128 (summand0/exp(1)) + exp (float128 (lsummand - 1))
                    else:
                        if not summand1:
                            summand1 = float128 (lsummand)
                        else:
                            summand1 = float128 (summand1 + log (1 + exp (float128 (lsummand - summand1))))
        if summand0 > 0:
            return -log (summand0) - 4500.0 * log (10) - divisor
        return -summand1 - 4500.0 * log (10)

def main():
    """
    main function
    """
    infile = 'bci.txt'
    abd = Abundance (infile)
    x = [48.1838, 1088.19]
    lnl = abd.etienne_likelihood (x)
    print lnl

if __name__ == "__main__":
    exit(main())
