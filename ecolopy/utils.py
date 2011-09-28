#!/usr/bin/python
"""
13 Jul 2011

some utils for abundance, essentially for computation of K(D,A)
"""

__author__  = "Francois-Jose Serra"
__email__   = "francois@barrabin.org"
__licence__ = "GPLv3"
__version__ = "0.12"

from gmpy2 import log, mul, mpfr, div, lngamma
from sys import stdout

global STIRLINGS
STIRLINGS = {}


def generate_random_neutral_abundance (model_name, size, **kwargs):
    '''
    :arguments model_name: model name (ewens, etienne, lognorm)
    
    :arguments size: size of the community (J)

    :returns: random neutral distribution of abundance
    
    other args should be of kind theta, I, m
    eg:
      import ecolopy
      ecolopy.generate_random_neutral_abundance('ewens', 100, theta=12, I=12)

    '''
    # import inside because otherwise strange never-ending import...
    from ecolopy.ecological_model import EcologicalModel
    model = EcologicalModel (model_name, **kwargs)
    return model.rand_neutral(size)


def lpoch (z, m):
    '''
    returns log Pochhammer symbol taking advantage of:
    Pochhammer symbol (z)_m = (z)(z+1)....(z+m-1) = gamma(z+m)/gamma(z)
    '''
    return lngamma(z+m) - lngamma(z)


def shannon_entropy(abund, inds):
    '''
    computes Shannon entropy (H) for a given abundance table
    and number of individues.

    :arguments abund: distribution of abundances as list

    :returns: shannon entropy
    
    '''
    shannon = mpfr(0.)
    for spe in abund:
        shannon -= spe * log (spe)
    return (shannon + inds * log (inds)) / inds


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
    return [counts[x] for x in sorted (counts)]


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
    to test multiplication of pylnomes try equality of the two functions:
            _mul_uneq_polyn(polyn_a, polyn_b, len_a, len_b)
                                ==
            _mul_simil_polyn(polyn_a, polyn_b, len_a, len_b)
    '''
    if not polyn_a:
        return polyn_b
    if not polyn_b:
        return polyn_a
    # set default values
    len_a = len (polyn_a)
    len_b = len (polyn_b)
    diff = abs (len_a - len_b)
    if len_a >= len_b:
        polyn_b = polyn_b + [mpfr(0)] * (diff)
    else:
        _  = polyn_a + [mpfr(0.)] * (diff)
        polyn_a = polyn_b[:]
        polyn_b = _
        len_a = len (polyn_a)
        len_b = len (polyn_b)
    # switcher
    if len_a > len_b*2: # most
        return _mul_uneq_polyn(polyn_a, polyn_b, len_a, len_b)
    return _mul_simil_polyn(polyn_a, polyn_b, len_a, len_b)


def _mul_simil_polyn(polyn_a, polyn_b, len_a, len_b):
    '''
    fast polynomial multiplication when polynomes are nearly equal
    -> iterates over factors
    '''
    def mult_one (la, lb, stop, start=0):
        '''
        internal that computes row multiplication of 2 lists
        start and stops are here to skip multiplications by zero
        '''
        return [mul (la[i], lb[i]) for i in xrange (start, stop)]
    max_len = len_a + len_b
    diff = len_a - len_b
    new = []
    for i in xrange (1, len_a +1):
        new.append (sum (mult_one (polyn_a[:i], polyn_b[i-1::-1], i)))
    len_a2 = len_a * 2 - 1
    len_a1 = len_a + 1
    len_a3 = len_a - 1
    for i in xrange (len_a, max_len - 1):
        new.append (sum (mult_one (polyn_a[i-len_a3 : len_a1],
                              polyn_b[len_a    : i-len_a :-1], len_a2-i, diff)))
    return new    


def _mul_uneq_polyn(polyn_a, polyn_b, len_a, len_b):
    '''
    fast polynomial multiplication when 1 polynome >~ 2 times larger.
    -> iterates over coefficients
    '''
    new = [mpfr(0)] * (len_a + len_b - 1)
    for i in xrange (len_a):
        pai = polyn_a[i]
        for j in xrange (len_b):
            new [i + j] += pai * polyn_b[j]
    return new


def pre_get_stirlings(max_nm, needed, verbose=True):
    '''
    takes advantage of recurrence function:
      s(n, m) = s(n-1, m-1) - (n-1) * s(n-1, m)
    and as  s(whatever, 0) = 0 :
      s(n+1, 1) = -n * s(n, 1)
    keep only needed stirling numbers (necessary for large communities)
    '''
    STIRLINGS [1, 1] = mpfr (1)
    for one in xrange (1, max_nm):
        if verbose and not one % 1000:
            stdout.write('\r    %s of %s, size: %s' % (one, max_nm,
                                                  STIRLINGS.__sizeof__()))
            stdout.flush()
        STIRLINGS [one+1, 1] = -one * STIRLINGS [one, 1]
        for two in xrange (2, one+1):
            STIRLINGS[one+1, two] = STIRLINGS[one, two-1] - one * STIRLINGS[one, two]
        STIRLINGS[one+1, one+1] = STIRLINGS[one, one]
        if one-1 in needed: continue
        for j in xrange (1, one):
            del (STIRLINGS[one-1,j])
    if verbose:
        stdout.write('\n')

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


