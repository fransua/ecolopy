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
from gmpy2 import log, mul, mpfr, gamma, div, lngamma

global STIRLINGS
STIRLINGS = {}

def poch (z, m):
    '''
    returns Pochhammer symbol taking advantage of:
    Pochhammer symbol (z)_m = (z)(z+1)....(z+m-1) = gamma(z+m)/gamma(z)
    '''
    return div (gamma(z+m), gamma(z))

def lpoch (z, m):
    '''
    returns log Pochhammer symbol taking advantage of:
    Pochhammer symbol (z)_m = (z)(z+1)....(z+m-1) = gamma(z+m)/gamma(z)
    '''
    return lngamma(z+m) - lngamma(z)

def table (out, spp=None):
    '''
    data to contingency table
    any kind of data
    '''
    if spp == None:
        spp = int (max (out))
    counts = dict (zip (set (out), [mpfr(0.)]*spp))
    for ind in out:
        counts[ind] += mpfr(1)
    return counts.values()

def factorial_div (one, two):
    '''
    computes a!/b!
    '''
    if one < two:
        return div(1., reduce (mul, xrange(one, two)))
    elif one > two:
        return reduce (mul, xrange(two, one))
    else:
        return mpfr (1.)

def mul_polyn(polyn_a, polyn_b):
    '''
    returns product of 2 polynomes
    '''
    if not polyn_a:
        return polyn_b
    if not polyn_b:
        return polyn_a
    # set default values
    len_a = len (polyn_a)
    len_b = len (polyn_b)
    polyn2 = [mpfr(0)] * (len_a + len_b - 1)
    for i in xrange (len_a):
        pai = polyn_a[i]
        for j in xrange (len_b):
            polyn2 [i + j] += pai * polyn_b[j]
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

def get_kda (abund, verbose=True):
    '''
    compute kda according to etienne formula
    '''

    abund = sorted (abund)
    abund = [mpfr(x) for x in abund]
    specabund = [sorted (list (set (abund))), table (abund)]
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
        # get of polyn1 exponential the number of individues for current species
        polyn2 = polyn1[:]
        for _ in xrange (1, specabund[1][i]):
            polyn1 = mul_polyn (polyn1, polyn2)
        # polyn = polyn * polyn1
        polyn = mul_polyn (polyn, polyn1)
    kda = []
    for i in polyn:
        kda.append (log (i))
    return kda




