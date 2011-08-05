#!/usr/bin/python
"""
13 Jul 2011

some utils for abundance
"""

__author__  = "Francois-Jose Serra"
__email__   = "francois@barrabin.org"
__licence__ = "GPLv3"
__version__ = "0.0"

from nzmath.combinatorial import stirling1
from gmpy2 import log, mul, mpfr
from itertools import combinations, product


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


def factorial_div (a, b):
    '''
    computes a!/b!
    '''
    if a < b:
        return 1. / reduce (mul, xrange(a, b))
    elif a > b:
        return reduce (mul, xrange(b, a))
    else:
        return mpfr (1.)


def exp_polyn(polyn, expon):
    '''
    computes polynomial exponential
    '''
    if expon <= 1:
        return polyn
    polyn2 = {}
    for it in combinations(polyn.keys(), expon-1):
        polyn2.setdefault (sum (it), 0)
        polyn2[sum (it)] = polyn2[sum (it)] + reduce (mul, (polyn[i] for i in it))
    return polyn2


def mul_polyn(polyn_a, polyn_b):
    '''
    returns product of polynomes
    '''
    if not polyn_a:
        return polyn_b
    polyn2 = {}
    for it in product(polyn_a.keys(), polyn_b.keys()):
        polyn2.setdefault(it[0] + it[1], 0)
        polyn2[it[0] + it[1]] = polyn2[it[0] + it[1]] + mul (polyn_a [it[0]],
                                                             polyn_b[it[1]])
    return polyn2
    
        
def get_kda (abund):
    '''
    compute kda according to etienne formula
    '''
    abund.sort ()
    specabund = [sorted (list (set (abund))), table (abund)]
    Sdiff     = len (specabund [1])
    polyn     = {}
    for i in xrange (Sdiff):
        polyn1 = {}
        print "  Computing species %s out of %s" % (i+1, Sdiff)
        for k in xrange (1, specabund[0][i]+1):
            coeff = abs (stirling1 (specabund[0][i], k)) * \
                    factorial_div (k, specabund[0][i])
            polyn1[k-1] = coeff
        # get of polyn1 exponential the number of individues for current species
        polyn1 = exp_polyn(polyn1, specabund[1][i])
        # polyn = polyn * polyn1
        polyn = mul_polyn (polyn, polyn1)
    kda = []
    for i in sorted(polyn, reverse=True):
        kda.append (float (log (polyn[i])))
    return kda



#try:
#    from math import lgamma
#except ImportError:
#    from scipy.special import gammaln as lgamma
# from numpy import logaddexp
#
#def get_kda_old (abund):
#    '''
#    slow function to get K(D,A)
#    '''
#    unq_abd = sorted (list (set (abund))) # g
#    len_unq = len (unq_abd)                    # NDA = i 
#    frq_unq = table (abund)               # f
#    max_a   = max(abund)                  # MaxA
#    phi = [0] * (max_a +1)
#    for s in xrange (len (abund)):
#        phi [abund[s]] += 1
#    T = [[] for _ in xrange (len_unq)]
#    T[0] = [0] * (unq_abd[0] + 1)
#    T[0][0] = 0
#    T[0][1] = 1
#    if unq_abd [0] != 1:
#        ls2 = [0] * (unq_abd[0] + 1)
#        ls2[1] = 1
#        for n in xrange (2, unq_abd [0]+1):
#            ls1 = [0] * (n + 1)
#            for im in xrange (n):
#                ls1 [im] = ls2 [im]
#            ls1 [n] = 0
#            for im in xrange (2, n+1):
#                ls2[im] = ls1[im] + float (ls1[im-1] * (im - 1)) / (n - 1)
#        for im in xrange (2, unq_abd[0]+1):
#            T[0][im] = ls2[im]
#    for _in in xrange (1, len_unq):
#        T[_in] = [0] * (unq_abd [_in] + 1)
#        T[_in] [1] = 1
#        ls2 = [0] * (unq_abd[_in] + 1)
#        for im in xrange (unq_abd [_in-1]+1):
#            ls2[im] = T[_in-1][im]
#        for n in xrange (unq_abd[_in-1]+1, unq_abd[_in]+1):
#            ls1 = [0] * (n + 1)
#            for im in xrange (n):
#                ls1 [im] = ls2 [im]
#            ls1 [n] = 0
#            for im in xrange (2, n + 1):
#                ls2 [im] = ls1 [im] + float (ls1 [im - 1] * (im - 1)) / (n - 1)
#        for im in xrange (2, unq_abd [_in] + 1):
#            T [_in][im] = ls2 [im]
#    old_T = T[:]
#    for i in xrange (len (old_T)):
#        for j in xrange (len (old_T[i])):
#            T[i][j] = log (old_T[i][j]) if old_T[i][j]>0 else float ('-inf')
#    K = [None] * (sum (abund) + 1)
#    poly2 = K[:]
#    K [0] = 0.0
#    degree  = 1
#    const = log (10**mpfr(4500.0/len (abund)))
#    for i in xrange (len_unq):
#        local_t = T [i]
#        deg = unq_abd [i]
#        for j in xrange (frq_unq[i]):
#            nn = 0
#            while K[nn] is None:
#                nn += 1
#            nn_start = nn + 1
#            for mm in xrange (1, deg + 1):
#                poly2 [nn+mm] = local_t[mm] + K[nn]
#            for nn in xrange (nn_start, degree):
#                for mm in xrange (1, deg):
#                    poly2 [nn+mm] = logaddexp (poly2 [nn+mm], local_t[mm] + K[nn])
#                poly2 [nn+mm+1] = local_t[mm+1] + K[nn]
#            degree += deg
#            for nn in xrange (nn_start):
#                K [nn] = None
#            for nn in xrange (nn_start, degree):
#                K [nn] = poly2[nn] - const
#                poly2[nn] = None
#    return K [len (abund):]
