#!/usr/bin/python
"""
11 Aug 2011


"""

__author__  = "Francois-Jose Serra"
__email__   = "francois@barrabin.org"
__licence__ = "GPLv3"
__version__ = "0.0"

"""
This script tests fmin_slsqp using Example 14.4 from Numerical Methods for
Engineers by Steven Chapra and Raymond Canale.  This example maximizes the
function f(x) = 2*x*y + 2*x - x**2 - 2*y**2, which has a maximum at x=2,y=1.
"""

from gmpy2 import mpfr, log, exp, lngamma
from abundance import Abundance
from utils import lpoch

def testfunc(x, abd, verbose=False):
    """
    """
    divisor     = newdivisor = sum0 = sum1 = mpfr(0.0)
    logx1       = log (x[1])
    mpfr11300   = mpfr (11300)
    lnum        = exp (mpfr11300)
    minus_11333 = mpfr(-11333.2)
    exp1        = exp (1)
    mpfr1       = mpfr(1)
    x0          = x[0]
    if not abd.factor: # define it
        abd.ewens_likelihood()
    poch1 = abd.factor + log (x0) * abd.S - lpoch (x[1], abd.J) + logx1 * abd.S + lngamma(x0)
    # blablabla
    logx0S   = log(x0 + abd.S)
    lgam_x0S = lngamma (x0 + abd.S)
    for A in xrange (abd.J - abd.S):
        lsum = poch1 + abd.K [A] + A * logx1 - lgam_x0S - divisor
        lgam_x0S += logx0S
        if sum0 < mpfr11300:
            if lsum < minus_11333:
                sum1 += log (mpfr1 + exp (lsum - sum1)) if sum1 else lsum
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

abd = Abundance('bci_short.txt')
abd.load_kda ('kda_short.pik')

for i in xrange(100):
    testfunc([abd.theta, abd.m], abd)

print 'etienne log lnL  :', testfunc([abd.theta, abd.m], abd)
print 'should be arround: -8752.7501182700435\n'
