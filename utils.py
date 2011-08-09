#!/usr/bin/python
"""
13 Jul 2011

some utils for abundance, essentially for computation of K(D,A)
"""

__author__  = "Francois-Jose Serra"
__email__   = "francois@barrabin.org"
__licence__ = "GPLv3"
__version__ = "0.0"

from nzmath.combinatorial import stirling1
from gmpy2 import log, mul, mpfr
from itertools import product

global STIRLINGS
STIRLINGS = {}

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

def factorial_div (one, two):
    '''
    computes a!/b!
    '''
    if one < two:
        return 1. / reduce (mul, xrange(one, two))
    elif one > two:
        return reduce (mul, xrange(two, one))
    else:
        return mpfr (1.)

def mul_polyn(polyn_a, polyn_b):
    '''
    returns product of polynomes
    '''
    if not polyn_a:
        return polyn_b
    if not polyn_b:
        return polyn_a
    # set default values
    polyn2 = {}
    for i in xrange (min(polyn_a.keys()) + min(polyn_b.keys()),
                     max(polyn_a.keys()) + max(polyn_b.keys()) + 1):
        polyn2[i] = mpfr(0)
    # compute comulative product
    for one, two in product (polyn_a.keys(), polyn_b.keys()):
        polyn2 [one+two] += mul (polyn_a [one], polyn_b [two])
    return polyn2

def pre_get_stirlings(max_nm):
    '''
    takes advantage of recurrence function:
    s(n,m) = s(n-1, m-1) - (n-1) * s(n-1,m)
    '''
    for one in xrange (1, max_nm+1):
        STIRLINGS [one, 1] = stirling1 (one, 1)
    for one in xrange (2, max_nm+1):
        for two in xrange (2, max_nm+1):
            if two > one:
                continue
            if two == one:
                STIRLINGS[one, two] = STIRLINGS [one-1, two-1] - 0
            else:
                STIRLINGS[one, two] = STIRLINGS [one-1, two-1] - \
                                          mul ((one-1), STIRLINGS [one-1,
                                                                       two])

def stirling (one, two):
    '''
    returns log unsingned stirling number, taking advantage of the fact that
    if x+y if odd signed stirling will be negative.
    takes also advantage of recurrence function:
    s(n,m) = s(n-1, m-1) - (n-1) * s(n-1,m)
    '''
    # return abs of stirling number
    if (one + two)%2:
        return -STIRLINGS [one, two]
    return STIRLINGS [one, two]

def get_kda (abund, verbose=False):
    '''
    compute kda according to etienne formula
    '''
    abund.sort ()
    specabund = [sorted (list (set (abund))), table (abund)]
    sdiff     = len (specabund [1])
    polyn     = {}
    pre_get_stirlings (max (specabund[0]))
    for i in xrange (sdiff):
        polyn1 = {0: mpfr(1.)}
        if verbose:
            print "  Computing species %s out of %s" % (i+1, sdiff)
        for k in xrange (1, specabund[0][i] + 1):
            coeff = stirling (specabund[0][i], k) * \
                    factorial_div (k, specabund[0][i])
            polyn1[k-1] = coeff
        #polyn1[k] = mpfr (1.)
        # get of polyn1 exponential the number of individues for current species
        polyn2 = polyn1.copy()
        for _ in xrange (1, specabund[1][i]):
            polyn1 = mul_polyn(polyn1, polyn2)
        # polyn = polyn * polyn1
        polyn = mul_polyn (polyn, polyn1)
    kda = []
    for i in sorted(polyn, reverse=True):
        kda.append (float (log (polyn[i])))
    return kda
