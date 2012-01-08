#!/usr/bin/python
"""
04 Jan 2012


"""

__author__  = "Francois-Jose Serra"
__email__   = "francois@barrabin.org"
__licence__ = "GPLv3"
__version__ = "0.0"

from gmpy2 import mpfr, mul
from random import random
from time import time

pola = [mpfr (random()*10000000000) for _ in xrange(1000)]
polb = [mpfr (random()*10000000000) for _ in xrange(4000)]
lena = len (pola)
lenb = len (polb)
mpf0 = mpfr(0)

def test_1(pola, polb,lena, lenb):
    new = [mpfr(0)]*(lena+lenb-1)
    for i in xrange(lena):
        for j in xrange(lenb):
            new[i+j] += pola[i] * polb[j]
    return new

def test_3(pola, polb,lena, lenb):
    new = [mpfr(0)]*(lena+lenb-1)
    for i in xrange(lena):
        for j in xrange(lenb):
            new[i+j] += pola[i] * polb[j]
    return new

def main():
    t0 = time()
    one = test_1(pola, polb, lena, lenb)
    print time()-t0
    t0 = time()
    tre = test_3(pola, polb, lena, lenb)
    print time()-t0

    print one==tre
    print one[:20]
    print ''
    print tre[:20]
    
exit (main())
            
