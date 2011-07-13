#!/usr/bin/python
"""
13 Jul 2011

some utils for abundance
"""

__author__  = "Francois-Jose Serra"
__email__   = "francois@barrabin.org"
__licence__ = "GPLv3"
__version__ = "0.0"

try:
    from math import lgamma
except ImportError:
    from scipy.special import gammaln as lgamma
from numpy import log, logaddexp

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

def get_kda (abund):
    '''
    slow function to get K(D,A)
    '''
    unq_abd = sorted (list (set (abund))) # g
    len_unq = len (unq_abd)                    # NDA = i 
    frq_unq = table (abund)               # f
    max_a   = max(abund)                  # MaxA
    phi = [0] * (max_a +1)
    for s in xrange (len (abund)):
        phi [abund[s]] += 1
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

    K = [None] * (sum (abund) + 1)
    poly2 = K[:]
    K [0] = 0.0
    degree  = 1
    const = log (10**(4500.0/len (abund)))
    for i in xrange (len_unq):
        local_t = T [i]
        deg = unq_abd [i]
        for j in xrange (frq_unq[i]):
            nn = 0
            while K[nn] is None:
                nn += 1
            nn_start = nn + 1
            for mm in xrange (1, deg + 1):
                poly2 [nn+mm] = local_t[mm] + K[nn]
            for nn in xrange (nn_start, degree):
                for mm in xrange (1, deg):
                    poly2 [nn+mm] = logaddexp (poly2 [nn+mm], local_t[mm] + K[nn])
                poly2 [nn+mm+1] = local_t[mm+1] + K[nn]
            degree += deg
            for nn in xrange (nn_start):
                K [nn] = None
            for nn in xrange (nn_start, degree):
                K [nn] = poly2[nn] - const
                poly2[nn] = None
    return K [len (abund):]
