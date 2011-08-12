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

from scipy import optimize
from numpy import array
from gmpy2 import mpfr, log, exp, lngamma
from abundance import Abundance
from utils import lpoch



def testfunc(x, abd, verbose=False):
    """
    """
    divisor     = newdivisor = sum0 = mpfr(0.0)
    logx1       = log (x[1])
    mpfr11300   = mpfr (11300)
    lnum        = exp (mpfr11300)
    minus_11333 = mpfr(-11333.2)
    exp1        = exp (1)
    mpfr1       = mpfr(1)
    x0          = x[0]
    if not abd.factor: # define it
        abd.ewens_likelihood()
    poch1 = abd.factor + log (x0) * abd.S - lpoch (x[1], abd.J) + logx1 * abd.S
    lgam_x0 = lngamma(x0)
    # blablabla
    def get_lsum(A, D):
        return poch1 + abd.K [A] + A * logx1 - \
               (lngamma (x0 + A + abd.S) - lgam_x0) - D
    sum1 = get_lsum (0,0)
    for A in xrange (1, abd.J - abd.S):
        lsum = get_lsum (A, divisor)
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

from time import time

print '\n\n'

abd = Abundance('bci.txt')
abd.load_kda ('kda_short.pik')


bounds = [(1, abd.S), (0, abd.J_tot)]

print "fmin_sequential least square"
t0 = time()
x = optimize.fmin_slsqp(testfunc, [abd.theta, abd.m], args=(abd, ),
#                        bounds = bounds,
                        iprint=0, full_output=1)
print "Elapsed time:", (time()-t0), "s"
print "Results",x
print "\n\n"

t0 = time()
print "fmin_sequential l-bfgs-b"
x = optimize.fmin_l_bfgs_b (testfunc, [abd.theta, abd.m], args=(abd, ),
                            iprint=0,
#                            bounds=bounds,
                            approx_grad=True)
print "Elapsed time:", (time()-t0), "s"
print "Results",x
print "\n\n"



t0 = time()
print "fmin_sequential tnc"
x = optimize.fmin_tnc (testfunc, [abd.theta, abd.m], args=(abd, ),
#                       bounds=bounds,
                       approx_grad=True)
print "Elapsed time:", (time()-t0), "s"
print "Results",x
print "\n\n"


t0 = time()
x = optimize.fmin (testfunc, [abd.theta, abd.I], args=(abd, ))
print "Elapsed time:", (time()-t0), "s"
print "Results",x

