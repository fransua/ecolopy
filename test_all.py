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

def test_ewens (kind, infile):
    '''
    test estimation of theta by max likelihood with ewens formula
    '''
    t0 = time ()
    abd = Abundance (infile)
    print ' Testing Ewens algorithm with BCI %s dataset' % kind
    if kind == 'full':
        wanted_theta = 32.214
        wanted_lnl   = 282.14
        time_max = 0.05
    elif kind == 'short':
        wanted_theta = 33.302
        wanted_lnl   = 162.742
        time_max = 0.05
    print '  -> Optimal value of theta: %.3f' % abd.theta
    if round (abd.theta, 3) != wanted_theta:
        exit ('\n test failed in ewens test (theta should have been %s)\n' %\
              wanted_theta)

    print '  -> likelihood of theta: %.3f' % abd.ewens_likelihood()
    if round (abd.ewens_lnl, 3) != wanted_lnl:
        exit ('\n test failed in ewens test (theta should have been %s)\n' %\
              wanted_lnl)
   
    print '\n  Elapsed time (should be < %s): %s sec\n' % (time_max, time() - t0)
    return abd

def test_etienne (kind, abd):
    '''
    test etienne basic functions
    TODO: replace values of x by maximum likelihood estimations
    '''
    if kind == 'short':
        etienne_lnl = -10217.052
        time_max = '10 sec'
    elif kind == 'full':
        etienne_lnl = -10087.610
        time_max = '100 sec'
    print ' Testing Etienne algorithm with BCI %s dataset' % kind
    t0 = time()
    x = [78.9545, 161.016]
    lnl = abd._etienne_likelihood (x, verbose=False)
    print '  -> Etienne lnL computed : %.3f' % lnl
    if round (lnl, 3) != etienne_lnl:
        exit('\n test failed in etienne test\n')
    print '\n  Elapsed time (should be < %s): %s sec\n' % (time_max, time() - t0)


def main():
    """
    main function
    """
    kind = argv[1]
    
    kinds = ['short', 'full']
    if not kind in kinds:
        exit('\nERROR: kind should be one of:%s\n' % \
             ((' "%s"' * len (kinds)) % tuple (kinds)))

    if kind == 'short':
        # reload abundance with shorter dataset
        infile = 'bci_short.txt'
    elif kind == 'full':
        infile = 'bci.txt'
    else:
        exit()

    print '\nstarting tests...\n\n'
    print '************************************************************\n'
    abd = test_ewens (kind, infile)
    print '************************************************************\n'
    test_etienne (kind, abd)
    print '************************************************************\n'


if __name__ == "__main__":
    exit(main())
