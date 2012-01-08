#!/usr/bin/python
"""
04 Jan 2012


"""

__author__  = "Francois-Jose Serra"
__email__   = "francois@barrabin.org"
__licence__ = "GPLv3"
__version__ = "0.0"


from gmpy2 import mpfr
from random import random

a = [mpfr(random()*100000) for _ in xrange(10)]

b = 'hola mundo'


