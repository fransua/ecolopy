'''
test efficiency of polynome multiplication
'''

from sys import argv
from random import random

def mul_polyn (la,lb):
    '''
    fast polynomial multiplication when polynomes are nearly equal
    '''
    def mult_one (la, lb, stop, start=0):
        return [la[i] * lb[i] for i in xrange (start, stop)]
    len_a = len (la)
    len_b = len (lb)
    diff = abs (len_a - len_b)
    if len_a >= len_b:
        lb = lb + [0.] * (diff)
    else:
        lc = la + [0.] * (diff)
        la = lb[:]
        lb = lc
    max_len = len_a + len_b
    new = []
    for i in xrange (1, len_a +1):
        new.append (sum (mult_one (la[:i], lb[i-1::-1], i)))
    len_a2 = len_a * 2 - 1
    len_a1 = len_a + 1
    len_a3 = len_a - 1
    for i in xrange (len_a, max_len - 1):
        new.append (sum (mult_one (la[i-len_a3 : len_a1],
                              lb[len_a    : i-len_a :-1], len_a2-i, diff)))
    return new

def rec_mult2 (la,lb):
    '''
    fast polynomial multiplication when 1 polynome >~ 2 times larger.
    '''
    len_a = len (la)
    len_b = len (lb)
    new = [0.] * (len_a + len_b - 1)
    for i in xrange (len_a):
        pai = la[i]
        for j in xrange (len_b):
            new[i+j] += pai*lb[j]
    return new

def main():

    test = argv[1]

    la = [random() for _ in xrange (5002)]
    lb = [random() for _ in xrange (3002)]

    if test == '1':
        c = mul_polyn(la,lb)
    elif test == '2':
        d = rec_mult2(la,lb)

    # something to quickly test that both function are returning the same:
    la = [random() for _ in xrange (81)]
    lb = [random() for _ in xrange (21)]
    c = mul_polyn(la,lb)
    d = rec_mult2(la,lb)
    if c == d:
        print '\n     test ok'
    else:
        print '\n     test NOT ok'
        print c
        print d


if __name__ == "__main__":
    exit(main())


