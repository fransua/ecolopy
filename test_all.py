#!/usr/bin/python
"""
12 Aug 2011

test functionality of abundance class
"""

__author__  = "Francois-Jose Serra"
__email__   = "francois@barrabin.org"
__licence__ = "GPLv3"
__version__ = "0.0"

from abundance import Abundance
from time import time
from sys import argv

def test_etienne (kind):
    '''
    test etienne basic functions
    TODO: replace values of x by maximum likelihood estimations
    '''
    infile = 'bci.txt'
    abd = Abundance (infile)
    if kind == 'short':
        etienne_lnl = '-10217.051548359694'
        time_max = '15 sec'
        abd = Abundance(abd.abund[:162])
    elif kind == 'full':
        etienne_lnl = '-10087.610335418967'
        time_max = '100 sec'
    print '\n Testing Etienne algorithm with BCI %s dataset\n\n' % kind
    t0 = time()
    x = [78.9545, 161.016]
    lnl = abd.etienne_likelihood (x, verbose=True)
    print ' Etienne lnL computed :', lnl
    print ' Etienne lnL should be:', etienne_lnl
    if str (lnl) != etienne_lnl:
        exit('\n test failed in etienne test\n')
    print '\n test ok\n'
    print ' Elapsed time (should be < %s): %s sec\n' % (time_max, time() - t0)


def main():
    """
    main function
    """
    kind = argv[1]
    kinds = ['short', 'full']
    if not kind in kinds:
        exit('\nERROR: kind should be one of:%s\n' % \
             ((' "%s"' * len (kinds)) % tuple (kinds)))
    test_etienne(kind)


if __name__ == "__main__":
    exit(main())
