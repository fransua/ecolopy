#!/usr/bin/python
"""
04 Jul 2011


"""

__author__  = "Francois-Jose Serra"
__email__   = "francois@barrabin.org"
__licence__ = "GPLv3"
__version__ = "0.0"

from random        import random
from scipy.stats   import chisqprob
from scipy         import optimize
from gmpy2         import mul, mpfr, log, exp, div, lngamma
from cPickle       import dump, load

from utils         import *

class Abundance (object):
    '''
    Compute theta m I ect... for a given data sample
     * J_tot is the size of the metacommunity, if not defined = J * 3
    '''
    def __init__ (self, data, J_tot=None):
        if type (data) != list:
            self.data_path = data
            data = self.parse_infile ()
        self.abund     = [mpfr(x) for x in sorted (data [:])]
        self.J         = mpfr(sum (data))
        self.S         = mpfr(len (data))
        self.theta     = mpfr(self._ewens_optimal_theta ())
        self.m         = self.theta / self.J / mpfr(2)
        self.I         = self.m * mpfr (self.J - 1) / (1 - self.m)
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
        factor = lngamma (self.J + 1)
        phi = table (self.abund)
        phi += [0] * int (max (self.abund) - len (phi))
        for spe in xrange (self.S):
            factor -= log (max (1, self.abund[spe]))
        for spe in xrange (max (self.abund)):
            factor -= lngamma (phi[spe] + 1)
        self.ewens_lnl = mpfr(lpoch (self.theta,
                               self.J)) - log (self.theta) * \
                               self.S - factor
        self.factor    = factor
        return self.ewens_lnl

    def _ewens_optimal_theta (self):
        '''
        optimize theta using theta likelihood function
        '''
        theta_like = lambda x: -self._ewens_theta_likelihood (x)
        return optimize.golden(theta_like, brack=[.01/self.J, self.J])

    def _etienne_optimal_params (self):
        '''
        optimize theta and I using etienne likelihood function
        using scipy package, values that are closest to the one proposed
        by tetame, are raised by fmin function:
          * fmin    : theta = 48.1838; I = 1088.19 ; lnL = -10091.166; t = 316s
        but lower likelihoods are found using other algorithms:
         -> bounds = [(1, self.S), (0, self.J_tot)]
          * slsqp   : theta = 143.43 ; I = 82.41   ; lnL = -10088.943; t = 227s
          * l_bfgs_b: theta = 140.11 ; I = 84.37   ; lnL = -10088.941; t = 343s
          * tnc     : theta = 52.62  ; I = 729.34  ; lnL = -10090.912; t = 297s
         -> bounds = [(1, self.S), (0, 2 * self.J_tot)]
          * slsqp   : theta = 129.92 ; I = 92.18   ; lnL = -10088.902; t = 107s
          * l_bfgs_b: theta = 133.88 ; I = 89.47   ; lnL = -10088.919; t = 246s
          * tnc     : theta = 34.51  ; I = 10241.17; lnL = -10085.534; t = 429s
         -> no bounds
          * slsqp   : theta = 82.46  ; I = 190.73  ; lnL = -10088.589; t = 86s
          * l_bfgs_b: theta = 144.87 ; I = 82.02   ; lnL = -10088.942; t = 220s
          * tnc     : theta = 688.55 ; I = 38.43   ; lnL = -10083.546; t = 429s
        '''
        bounds = [(1, self.S), (0, self.J_tot)]
        args_like = lambda x, y: -self._etienne_likelihood (x, y)
        return optimize.fmin (args_like, [self.theta, self.I], args=(False, ))

    def _ewens_theta_likelihood (self, theta):
        '''
        returns the likelihood of theta for a given dataset
        '''
        if theta < 0:
            return mpfr ('-inf')
        return self.S * log (theta) + lngamma (theta) - lngamma (theta + self.J)

    def _shannon_entropy (self):
        '''
        computes shannon entropy
        '''
        H = mpfr(0.)
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
            H = mpfr(0.)
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

    def _etienne_likelihood(self, x, verbose=True):
        '''
        logikelihood function where
          x[0] = theta
          x[1] = I
          K[A] = log (K(D,A))
          K[A] = K(D,A) in Etienne paper
        cf Etienne, 2005
        x = [48.1838, 1088.19]
        returns log likelihood of given theta and I
        '''
        divisor     = newdivisor = sum0 = mpfr(0.0)
        logx1       = log (x[1])
        mpfr11300   = mpfr (11300)
        lnum        = exp (mpfr11300)
        minus_11333 = mpfr(-11333.2)
        exp1        = exp (1)
        mpfr1       = mpfr(1)
        x0          = x[0]
        if not self.factor: # define it
            self.ewens_likelihood()
        lgam_x0 = lngamma(x0)
        # blablabla
        if not self.factor: # define it
            self.ewens_likelihood()
        if not self.K:
            if verbose:
                print "\nGetting K(D,A) according to Etienne 2005 formula:"
            self._get_kda (verbose=verbose)
        # pre calculate som staff for get_lsum
        poch1 = self.factor + log (x0) * self.S - \
                lpoch (x[1], self.J) + logx1 * self.S
        def __get_lsum(A, D):
            return poch1 + self.K [A] + A * logx1 - \
                   (lngamma (x0 + A + self.S) - lgam_x0) - D
        sum1 = __get_lsum (0,0)
        for A in xrange (1, self.J - self.S):
            lsum = __get_lsum (A, divisor)
            if sum0 < mpfr11300:
                if lsum < minus_11333:
                    sum1 += log (mpfr1 + exp (lsum - sum1))
                else:
                    sum0 += exp (lsum)
                continue
            if lsum > mpfr11300:
                newdivisor = lsum - mpfr11300
                divisor   += newdivisor
                sum0       = sum0 / exp (newdivisor) + lnum
                continue
            divisor  += mpfr1
            sum0 = sum0 / exp1 + exp (lsum - mpfr1)
        if sum0 > 0:
            return -log (sum0) - 4500.0 * log (10) - divisor
        return -sum1 - 4500.0 * log (10)

    def dump_kda (self, outfile):
        '''
        save kda with pickle
        '''
        if not self.K:
            self._get_kda (verbose=False)
        dump (self.K, open (outfile, 'w'))

    def load_kda (self, infile):
        '''
        load kda with pickle from infile
        '''
        self.K = load (open (infile))

    def _get_kda (self, verbose=True):
        '''
        compute kda according to etienne formula
        '''
        specabund = [sorted (list (set (self.abund))), table (self.abund)]
        sdiff     = len (specabund [1])
        polyn     = []
        # compute all stirling numbers taking advantage of recurrence function
        pre_get_stirlings (max (specabund[0]))
        for i in xrange (sdiff):
            if verbose:
                print "  Computing species %s out of %s" % (i+1, sdiff)
            polyn1 = []
            for k in xrange (1, specabund[0][i] + 1):
                coeff = stirling (specabund[0][i], k) * \
                        factorial_div (k, specabund[0][i])
                polyn1.append (coeff)
            if not polyn1:
                polyn1.append(mpfr(1.))
            # polyn1 exponential the number of individues for current species
            polyn2 = polyn1[:]
            for _ in xrange (1, specabund[1][i]):
                polyn1 = mul_polyn (polyn1, polyn2)
            # polyn = polyn * polyn1
            polyn = mul_polyn (polyn, polyn1)
        kda = []
        for i in polyn:
            kda.append (log (i))
        self.K = kda
