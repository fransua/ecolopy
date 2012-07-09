#!/usr/bin/python
"""
22 Dec 2011


"""

__author__  = "Francois-Jose Serra"
__email__   = "francois@barrabin.org"
__licence__ = "GPLv3"
__version__ = "0.0"

try:
    from gmpy2 import mpfr, mul
except ImportError:
    warn("WARNING: GMPY2 library not found, using numpy")
    from numpy import float128 as mpfr
    from operator import mul

from sys import setrecursionlimit

setrecursionlimit(10000)

mpf0 = mpfr(0)

class Polynomial(object):
    """
    """
    def __init__(self, plist=None):
        """
        """
        if plist is None:
            plist = []
        elif plist == []:
            plist = [mpfr(1.0)]
        self.plist = [mpfr (p) for p in plist]

    def add(self, something):
        """
        similar to append for lists
        """
        self.plist.append (something)

    def __str__(self):
        string = []
        for i in xrange (len (self)):
            if not self.plist[i]: continue
            if i > 1:
                string.append ('{0}x^{1}'.format(int (self.plist[i]), i))
            elif i==1:
                string.append ('{0}x'.format(int (self.plist[i])))
            else:
                string.append ('{0}'.format(int (self.plist[i])))
        return ' + '.join (reversed (string))

    def __repr__(self):
        return str(self.plist)

    def __neg__(self):
        return self * Polynomial([-1])

    def __getitem__(self, key): 
        return self.plist[key]

    def __len__(self):
        return len (self.plist)
    
    def __add__(self, poly):
        return add (self, poly)

    def __iadd__ (self, poly):
        return self + poly

    def __sub__(self, poly):
        return self + (-poly)
        
    def __mul__ (self, poly):   
        return Polynomial (mul_polyn (self.plist, poly.plist))

    def __pow__ (self, p):
        #return Polynomial (reduce (mul_polyn, (self.plist for _ in xrange(p))))
        return poly_2_pow (self, int(p))

    def __iter__(self):
        return self.iter_elts()

    def iter_elts(self):
        for i in self.plist:
            yield i

def add(poly_a, poly_b):
    if not type(poly_b) is Polynomial:
        raise Exception('can only sum Polynomials')
    if len (poly_a) < len (poly_b):
        return Polynomial (poly_b) + poly_a
    new = poly_a.plist[:]
    for i in xrange (len (poly_b)):
        new[i] += poly_b[i]
    return Polynomial (new)


def poly_2_pow (poly, p, r=()):
    if p == 1:
        for i in r:
            poly = poly * i
        return poly
    if p == 0:
        return  Polynomial([0])
    if p%2:
        r = r + (Polynomial(poly.plist[:]),)
    poly = poly * poly
    p = p/2
    if p == 2:
        poly = poly * poly
        for i in r:
            poly = poly * i
        return poly
    return poly_2_pow (poly,p,r)


def mul_polyn(polyn_a, polyn_b):
    '''
    computes the product of 2 polynomials, depending of the differences in length
    of the two polynomials, this function will call one of:
    * _mul_uneq_polyn: when length of polyn_a is >= length of polyn_b, will iterate over coefficient.
    * _mul_simil_polyn: in case both polynomials have equal length, will iterate over factors.

    to test multiplication of pylnomials try equality of the two functions:
    mul_uneq_polyn(polyn_a, polyn_b, len_a, len_b) == _mul_simil_polyn(polyn_a, polyn_b, len_a, len_b)

    **Example:**
    ::

      from ecolopy_dev.utils import mul_polyn
      # (2 + 3^x + 5x^2)  * (x)
      mul_polyn([2,3,5], [0,1])
      # will return: [mpfr('0.0'), mpfr('2.0'), mpfr('3.0'), mpfr('5.0')]
      # that is: 0 + 2x + 3x^2 + 5x^3
      
    :argument polyn_a: list of indices of polynomial
    :argument polyn_b: list of indices of polynomial (e.g.: [1,3,5,0,2] for :math:`2 + 3^x + 5x^2 + 0x^3 + 2x^4`)
    :returns: a list representing multiplication of the two polynomials

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
        _  = Polynomial(polyn_a + [mpfr(0.)] * (diff))
        polyn_a = Polynomial(polyn_b[:])
        polyn_b = _
        len_a = len_b
        len_b = len (polyn_b)
    # switcher
    if len_a > len_b*2 and True: # most
        return _mul_uneq_polyn(polyn_a, polyn_b, len_a, len_b)
    return _mul_simil_polyn(polyn_a, polyn_b, len_a, len_b)


def _mul_simil_polyn(polyn_a, polyn_b, len_a, len_b):
    '''
    fast polynomial multiplication when polynomials are nearly equal
    -> iterates over factors
    '''
    def mult_one (la, lb, stop, start=0):
        '''
        internal that computes row multiplication of 2 lists
        start and stops are here to skip multiplications by zero
        '''
        return [mul(la[i], lb[i]) for i in xrange (start, stop)]
    max_len = len_a + len_b
    diff = len_a - len_b
    new = []
    for i in xrange (1, len_a +1):
        new.append (sum (mult_one (polyn_a[:i], polyn_b[i-1::-1], i)))
    len_a2 = len_a * 2 - 1
    len_a1 = len_a + 1
    len_a3 = len_a - 1
    for i in xrange (len_a, max_len - 1):
        new.append (sum (mult_one (polyn_a[i-len_a3:len_a1],
                                   polyn_b[len_a:i-len_a:-1], len_a2-i, diff)))
    return new


def _mul_uneq_polyn(polyn_a, polyn_b, len_a, len_b):
    '''
    fast polynomial multiplication when 1 polynomial >~ 2 times larger.
    -> iterates over coefficients
    '''
    new = [mpf0] * (len_a + len_b - 1)
    for ij, pa, pb in ((i+j, polyn_a[i], polyn_b[j]) for i in xrange(len_a) for j in xrange(len_b)):
        new [ij] += mul(pa,pb)
    #for i in xrange (len_a):
    #    pai = polyn_a[i]
    #    for j in xrange (len_b):
    #        new [i + j] += pai * polyn_b[j]
    return new
