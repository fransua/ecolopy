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

from gmpy2 import mpfr, log, exp, lngamma, gamma
from abundance import Abundance
from utils import lpoch
from sys import argv

def testfunc1 (x, abd, verbose=False):
    """
    """
    lik = mpfr(0.0)
    x0          = x[0] + abd.S
    I = float(x[1])/(1 - x[1]) * (abd.J - 1)
    logx1       = log (I)
    if not abd.factor: # define it
        abd.ewens_likelihood()
    poch1 = exp (abd.factor + log (x[0]) * abd.S - lpoch (I, abd.J) + \
                 logx1 * abd.S + lngamma(x[0]))
    gam_x0S = gamma (x0)
    for A in xrange (abd.J - abd.S):
        lik += poch1 * exp (abd.K [A] + A * logx1) / gam_x0S
        gam_x0S *= x0 + A
    return -log (lik)

abd = Abundance('bci.txt')
abd.load_kda ('kda.pik')

if argv[1] == '1':
    testfunc = testfunc1
else:
    testfunc = testfunc2
    
for theta in xrange (1, 100, 10):
    for I in xrange (1, 1000, 93):
        testfunc([theta, I], abd)

print 'etienne log lnL  :', testfunc1([abd.theta, abd.m], abd)
print 'etienne log lnL  :', testfunc2([abd.theta, abd.m], abd)
print 'etienne log lnL  :', testfunc1([abd.theta, abd.I], abd)
print 'etienne log lnL  :', testfunc2([abd.theta, abd.I], abd)
print 'should be arround: -8752.7501182700435\n'
