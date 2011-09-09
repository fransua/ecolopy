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
from sys import argv, stderr
from numpy import mean
from random_neutral import *
from utils import *

def test_ewens (kind, abd):
    '''
    test estimation of theta by max likelihood with ewens formula
    '''
    t0 = time ()
    print ' Testing Ewens algorithm with BCI %s dataset' % kind
    if kind == 'full':
        wanted_theta = 34.962
        wanted_lnl   = 318.849
        time_max = 0.05
    elif kind == 'short':
        wanted_theta = 33.302
        wanted_lnl   = 162.742
        time_max = 0.05
    abd.ewens_optimal_params()
    print '  -> Optimal value of theta: %.3f' % abd.params['ewens']['theta']
    if round (abd.params ['ewens']['theta'], 3) != wanted_theta:
        stderr.write ('\n test failed in ewens test (theta should have been %s)\n' %\
              wanted_theta)

    print '  -> likelihood of theta: %.3f' % abd.params ['ewens']['lnL']
    if round (abd.params ['ewens']['lnL'], 3) != wanted_lnl:
        stderr.write ('\n test failed in ewens test (theta should have been %s)\n' %\
              wanted_lnl)
   
    print '\n  Elapsed time (should be < %s): %s sec\n' % (time_max, time() - t0)

def test_etienne (kind, abd):
    '''
    test etienne basic functions
    '''
    if kind == 'short':
        etienne_lnl  = 144.581
        wanted_theta = 78.954
        wanted_m     = 0.03677
        time_max = '15 sec'
    elif kind == 'full':
        etienne_lnl = 308.725
        wanted_theta = 47.674
        wanted_m     = 0.09342
        time_max = '150 sec'
    print ' Testing Etienne algorithm with BCI %s dataset' % kind
    t0 = time()
    for test in ['fmin', 'slsqp', 'l_bfgs_b', 'tnc']:
        print test
        try:
            abd.etienne_optimal_params (method=test)
        except Exception as err:
            print err
            continue
        print '  -> Optimal value of theta: %.3f' % abd.params['etienne']['theta']
        if round (abd.params ['etienne']['theta'], 3) != wanted_theta:
            stderr.write ('\n test failed in ewens test (theta should have been %s)\n' %\
                  wanted_theta)
        print '  -> Optimal value of m: %.5f' % abd.params['etienne']['m']
        if round (abd.params ['etienne']['m'], 5) != wanted_m:
            stderr.write ('\n test failed in etienne test (m should have been %s)\n' %\
                  wanted_m)
        print '  -> Etienne lnL computed : %.3f' % abd.params ['etienne']['lnL']
        if round (abd.params ['etienne']['lnL'], 3) != etienne_lnl:
            stderr.write ('\n test failed in etienne test\n')
        print '\n  Elapsed time (should be < %s): %s sec\n' % (time_max, time() - t0)


def test_random_neutral (abd):
    '''
    '''
    thetas = []
    immigs = []
    for _ in xrange (100):
        print _
        new = Abundance (abd.rand_neutral (5.0,1.5, model='lognorm'))
        new.etienne_optimal_params()
        thetas.append (new.params['etienne']['theta'])
        immigs.append (new.params['etienne']['I'])
    print mean (thetas), '50'
    print mean (immigs), '1000'


def main():
    """
    main function
    kind = 'short'
    infile = 'bci_short.txt'
    """
    try:
        kind = argv[1]
    except IndexError:
        kind = None
    
    kinds = ['short', 'full']
    if not kind in kinds or not kind:
        exit('\nERROR: kind should be one of:%s\n' % \
             ((' "%s"' * len (kinds)) % tuple (kinds)))

    if kind == 'short':
        # reload abundance with shorter dataset
        abd = Abundance ('bci_short.txt')
    elif kind == 'full':
        abd = Abundance ('bci_full.txt')
    else:
        exit()

    print '\nstarting tests...\nwith %s BCI dataset:' % kind
    print abd
    print '************************************************************\n'
    print '  Dataset with J: %s, S: %s, H: %s\n' % (abd.J, abd.S, abd.shannon)
    print '************************************************************\n'
    test_ewens (kind, abd)
    print '************************************************************\n'
    test_etienne (kind, abd)
    print '************************************************************\n'
    print ' Testing load/dump data...',
    abd.dump_params('test_%s.pik' % kind)
    abd.load_params('test_%s.pik' % kind)
    # and again to test update
    abd.dump_params('test_%s.pik' % kind)
    print 'ok\n'
    print abd.params['ewens']['theta'], abd.params['ewens']['I']
    print abd.params['etienne']['theta'], abd.params['etienne']['I']
    print '************************************************************\n'
    print 'LRT between Ewens and Etienne model (1 df): ',
    print abd.lrt('ewens', 'etienne'), '\n'
    print '************************************************************\n'
    t0 = time()
    gens = 500
    print 'Neutrality test p-value (under Ewens model):',
    print abd.test_neutrality(model='ewens', gens=gens)
    print 'should be arround: %s' % (0.9 if kind=='full' else 1.0)
    print '%s generations computed in %ss' % (gens, time()-t0)
    t0 = time()
    print 'Neutrality test p-value (under Etienne model):',
    print abd.test_neutrality(model='etienne', gens=gens)
    print 'Lognormality test p-value (under Etienne model):',
    abd.lognorm_optimal_params()
    print abd.test_neutrality(model='lognorm', gens=gens)
    print 'lognorm lnL', abd.lognorm_likelihood()
    print 'ewens   lnL',abd.params['ewens']['lnL']
    print 'etienne lnL',abd.params['etienne']['lnL']
    print abd.params['factor']
    print lpoch (abd.params['etienne']['theta'], abd.J)
    print abd.params['factor'] - lpoch (abd.params['ewens']['theta'], abd.J) + abd.params['ewens']['lnL']
    print abd.params['factor'] - lpoch (abd.params['etienne']['theta'], abd.J) + abd.params['etienne']['lnL']
    print 'should be arround: %s' % (0.1 if kind=='full' else 0.01)
    print '%s generations computed in %ss' % (gens, time()-t0)
    print '************************************************************\n'
    abd._kda = None
    test_random_neutral(abd)
    print '\n\nAll test OK!\n'


if __name__ == "__main__":
    exit(main())
