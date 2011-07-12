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
from numpy       import log
from random      import random
from scipy.stats import chisqprob
from scipy       import optimize
from numpy import logaddexp

def table (out, spp=None):
    '''
    data to contingency table
    any kind of data
    '''
    if spp == None:
        spp = max (out)
    counts = dict (zip (set (out), [0]*spp))
    for ind in out:
        counts[ind] += 1
    return counts.values()

lpochham = lambda x, n: lgamma(x+n)-lgamma(x)

class Abundance (object):
    '''
    Compute theta m I ect... for a given data sample
    '''
    def __init__ (self, data, J_tot=None):
        if type (data) != list:
            self.data_path = data
            data = self.parse_infile ()
        self.abund = data [:]
        self.J     = sum (data)
        self.S     = len (data)
        self.theta = self._ewens_optimal_theta ()
        self.m     = self.theta / self.J / 2
        self.I     = self.m * float (self.J - 1) / (1 - self.m)
        self.J_tot = J_tot if J_tot else self.J * 3
        self.H     = self._shannon_entropy ()

    def _lrt_ewens_etienne(self):
        '''
        likelihood-ratio-test between two neutral models
        returns p-value of rejection of alternative model
        '''
        chisqprob(2 * (self.ewens_likelihood() - 0), df=1)

    def ewens_likelihood (self):
        '''
        get likelihood value of Ewens according to parthy/tetame
        returns likelihood and factor
        '''
        factor = lgamma (self.J+1)
        phi = table (self.abund)
        phi += [0] * (max (self.abund)-len (phi))
        for spe in xrange (self.S):
            factor -= log (max (1, self.abund[spe]))
        for spe in xrange (max (self.abund)):
            factor -= lgamma (phi[spe] + 1)
        return (lpochham (self.theta,
                          self.J) - log (self.theta) * self.S - factor,
                factor)

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
          x[0]=theta
          x[1]=I,K[A]=log(K(D,A))
          
        K[A]=K(D,A)in Etienne paper
        
        cf Etienne, 2005
        '''
        divisor = 0.0
        newdivisor = 0.0
        summand = 0.0
        for A in xrange (self.S, self.J):
            lsummand = factor + log(x[0]) * self.S - \
                       lpochham (x[1], int (self.J)) + K [int (A)] + \
                       A * log (x[1]) - lpochham (x[0], int (A)) - divisor

    def get_kda (self):

        unq_abd = sorted (list (set (self.abund))) # g
        len_unq = len (unq_abd)                    # NDA = i 
        frq_unq = table (self.abund)               # f
        max_a   = max(self.abund)                  # MaxA
        phi = [0] * (max_a +1)
        for s in xrange (self.S):
            phi [self.abund[s]] += 1
        T = [[] for _ in xrange (len_unq)]
        T[0] = [0] * (unq_abd[0] + 1)
        T[0][0] = 0
        T[0][1] = 1
        if unq_abd [0] != 1:
            ls2 = [0] * (unq_abd[0] + 1)
            ls2[1] = 1
            for n in xrange (2, unq_abd [0]+1):
                ls1 = [0] * (n + 1)
                for im in xrange (n):
                    ls1 [im] = ls2 [im]
                ls1 [n] = 0
                for im in xrange (2, n+1):
                    ls2[im] = ls1[im] + float (ls1[im-1] * (im - 1)) / (n - 1)
            for im in xrange (2, unq_abd[0]+1):
                T[0][im] = ls2[im]

        for _in in xrange (1, len_unq):
            T[_in] = [0] * (unq_abd [_in] + 1)
            T[_in] [1] = 1
            ls2 = [0] * (unq_abd[_in] + 1)
            for im in xrange (unq_abd [_in-1]+1):
                ls2[im] = T[_in-1][im]
            for n in xrange (unq_abd[_in-1]+1, unq_abd[_in]+1):
                ls1 = [0] * (n + 1)
                for im in xrange (n):
                    ls1 [im] = ls2 [im]
                ls1 [n] = 0
                for im in xrange (2, n + 1):
                    ls2 [im] = ls1 [im] + float (ls1 [im - 1] * (im - 1)) / (n - 1)
            for im in xrange (2, unq_abd [_in] + 1):
                T [_in][im] = ls2 [im]
        old_T = T[:]
        for i in xrange (len (old_T)):
            for j in xrange (len (old_T[i])):
                T[i][j] = log (old_T[i][j]) if old_T[i][j]>0 else float ('-inf')

        K = [None] * (self.J + 1)
        poly2 = K[:]
        K [0] = 0.0
        degree  = 1
        const = log (10**(4500.0/self.S))
        for i in xrange (len_unq):
            local_t = T [i]
            deg = unq_abd [i]
            for j in xrange (frq_unq[i]):
                nn = 0
                while K[nn] is None:
                    nn += 1
                for mm in xrange (1, deg + 1):
                    #print "-",nn, mm, i, j, degree
                    poly2 [nn+mm] = local_t[mm] + K[nn]
                for nn in xrange (nn+1, degree):
                    for mm in xrange (1, deg):
                        #print "*", nn, mm, i, j, degree
                        poly2 [nn+mm] = logaddexp (poly2 [nn+mm], local_t[mm] + K[nn])
                    #print "-",nn, mm, i, j, degree
                    poly2 [nn+mm+1] = local_t[mm+1] + K[nn]
                degree += deg
                nn = 0
                while poly2[nn] is None:
                    K [nn] = None
                    nn+=1
                for nn in xrange (nn, degree):
                    K [nn] = poly2[nn] - const
                    poly2[nn] = None
        K = K [self.S:]


#        K = [None] * (self.J + 1)
#        poly2 = K[:]
#        K [0] = 0.0
#        degree  = 1
#        const = log (10**(4500.0/self.S))
#        for i in xrange (len_unq):
#            for _ in xrange (frq_unq[i]):
#                for nn in xrange (degree):
#                    if K[nn] is None: continue
#                    for mm in xrange (1, unq_abd [i]+1):
#                        if poly2 [nn+mm] is None:
#                            poly2 [nn+mm] = T[i][mm] + K[nn]
#                        else:
#                            poly2 [nn+mm] = logaddexp (poly2 [nn+mm], T[i][mm] + K[nn])
#                degree += unq_abd[i]
#                for nn in xrange (degree):
#                    K [nn] = None if poly2[nn] is None else (poly2[nn] - const)
#                    poly2[nn] = None
#        K = K [self.S:]

def main():
    """
    main function
    """
    infile = 'bci.txt'
    abd = Abundance (infile)
    abd.get_kda()

if __name__ == "__main__":
    exit(main())

