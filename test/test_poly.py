#!/usr/bin/python
"""
02 Jan 2012


"""

from ecolopy_dev.polynomial import Polynomial, mul_polyn
from ecolopy_dev.utils import *
from gmpy2 import mpfr


poly = [1,0,2,3,6,5,8,65,3,2,5,2,3,4,5]*100

poly = [mpfr (p) for p in poly]

Poly = Polynomial(poly)

multipliers = ([2], [3,5,0,9], [1], [0], [0,-1], [1,0,2,3,6,5,8,65,3,2,5,2,3,4,5]*100, [3,598934354534643549434354943138499,6,98,7,4,89,8,9,8,9,8,9,8,9,8,6,5,5,3,4,6,8,4,3,5,256,3,3569874,256,3,12589,999,9,9,9,9,9,9,9,9,9,9,9,9,9,9,9,9,9,9,9,9,9,9,9,9,9,9,9,9,2,2,5,4,5,6,6,5,6,5,4,5,2,8,4,8,4,5,1,3,2,1,3,8,5,4,3,5,1,3,8,9,9,9,9,9,9,9,9,9,9,9,9,9,9,9,9,9,9,9,9,9,9,9,99999999999999999,9,9,9,9,9,9,9,2])
multipliers = ([999999999999999999999999999999999999999999999999999999999999999999999999,0,2,3,6,5,8,65,3,2,5,2,3,4,5]*100,)

# multiplication
for m in multipliers:
    m = [mpfr(i) for i in m]
    print m
    p = Poly * Polynomial(m)
    one = p.plist
    two = mul_polyn(poly, m)
    while two[-1]==0.0 and len (two)>1:
        two.pop(-1)
    while one[-1]==0.0 and len (one)>1:
        one.pop(-1)
    print Polynomial(one)
    print Polynomial(two)
    print one == two

poly = [1,0,2,3,6,5,8,65,3,2,5,2,3,4,5]
poly = [mpfr (p) for p in poly]
Poly = Polynomial(poly)
# powers

for pw in [1,2,5,8,10,20,24]:
    print 'power:', pw
    pw = mpfr(pw)
    p = Poly ** pw
    one = p.plist
    two = power_polyn(poly, pw)
    tre = Polynomial (reduce (mul_polyn, (Poly.plist for _ in xrange(pw)))).plist
    while tre[-1]==0.0 and len (tre)>1:
        tre.pop(-1)
    while two[-1]==0.0 and len (two)>1:
        two.pop(-1)
    while one[-1]==0.0 and len (one)>1:
        one.pop(-1)
    print Polynomial(one)
    print Polynomial(two)
    print Polynomial(tre)
    print one == two
    print one == tre
    print two == tre


