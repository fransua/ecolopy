#!/usr/bin/python
"""
12 Aug 2011

test functionality of abundance class
"""

__author__  = "Francois-Jose Serra"
__email__   = "francois@barrabin.org"
__licence__ = "GPLv3"
__version__ = "0.13"

from ecolopy_dev                import Community
from time                       import time
from sys                        import argv, stderr
from numpy                      import mean
from ecolopy_dev.models         import EcologicalModel
from ecolopy_dev.random_neutral import *
from ecolopy_dev.utils          import *


def test_ewens (kind, abd, test_failed):
    '''
    test estimation of theta by max likelihood with ewens formula
    '''
    t0 = time ()
    test_ok=True
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
    abd.set_current_model ('ewens')
    model = abd.get_model ('ewens')
    print '  -> Optimal value of theta: %.3f' % abd.theta
    if round (abd.theta, 3) != wanted_theta:
        stderr.write ('\n test failed in ewens test (theta should have been %s)\n' %\
              wanted_theta)
        test_ok=False

    print '  -> likelihood of theta: %.3f' % model.lnL
    if round (model.lnL, 3) != wanted_lnl:
        stderr.write ('\n test failed in ewens test (theta should have been %s)\n' %\
              wanted_lnl)
        test_ok=False
   
    print '\n  Elapsed time (should be < %s): %s sec\n' % (time_max, time() - t0)
    if not test_ok:
        test_failed.append ('etienne optimization')
    return test_failed


def test_etienne (kind, abd, test_failed):
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
    models = []
    test_ok = True
    t0 = time()
    #abd._get_kda()
    print '\n      Elapsed time (should be < %s): %s sec\n' % (time_max, time() - t0)
    print '  Testing different optimization strategies: '
    for test in ['fmin', 'slsqp', 'l_bfgs_b', 'tnc']:
        t0 = time()
        print '    ' + test
        try:
            abd.etienne_optimal_params (method=test)
        except Exception as err:
            print err
            continue
        model = abd.get_model('etienne')
        print '    -> Optimal value of theta: %.3f' % model.theta
        print '    -> Optimal value of m: %.5f' % model.m
        print '    -> Etienne lnL computed : %f' % model.lnL
        print '\n      Elapsed time (should be < %s): %s sec\n' % (time_max, time() - t0)
        models.append (model)
    model = min (models, key=lambda x: x.lnL)
    print '\nlikelihood of better optimization:', model.lnL
    print '  -> Optimal value of theta: %.3f' % model.theta
    if round (model.theta, 3) != wanted_theta:
        stderr.write ('\n test failed in ewens test (theta should have been %s)\n' %\
              wanted_theta)
        test_ok = False
    print '  -> Optimal value of m: %.5f' % model.m
    if round (model.m, 5) != wanted_m:
        stderr.write ('\n test failed in etienne test (m should have been %s)\n' %\
              wanted_m)
        test_ok = False
    print '  -> Etienne lnL computed : %f' % model.lnL
    if round (model.lnL, 3) != etienne_lnl:
        stderr.write ('\n test failed in etienne test\n')
        test_ok = False
    abd.set_model (min (models, key=lambda x: x.lnL))
    if not test_ok:
        test_failed.append ('etienne optimization')
    return test_failed


def test_random_neutral (abd):
    '''
    '''
    thetas = []
    immigs = []
    model = EcologicalModel ('lognorm', theta=0.5, I=1.5)
    for _ in xrange (100):
        print _
        new = Abundance (model.rand_neutral (abd.J))
        new.etienne_optimal_params()
        new.set_current_model ('etienne')
        thetas.append (new.theta)
        immigs.append (new.I)
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
        abd = Community ('../dataset_trial/bci_short.txt')
    elif kind == 'full':
        abd = Community ('../dataset_trial/bci_full.txt')
    else:
        exit()
    test_failed = []
    abd.fit_model(name='etienne')
    abd.fit_model(name='ewens')
    ## print '\nstarting tests...\nwith %s BCI dataset:' % kind
    ## print abd
    ## print '************************************************************\n'
    ## print '  Dataset with J: %s, S: %s, H: %s\n' % (abd.J, abd.S, abd.shannon)
    ## print '************************************************************\n'
    ## test_failed = test_ewens (kind, abd, test_failed)
    ## print '************************************************************\n'
    ## test_failed = test_etienne (kind, abd, test_failed)
    ## print '************************************************************\n'
    ## print ' Testing load/dump data...',
    ## abd.dump_abundance('test_%s.pik' % kind)
    ## abd.load_abundance('test_%s.pik' % kind)
    ## # and again to test update
    ## abd.dump_abundance('test_%s.pik' % kind)
    ## print 'ok\n'
    ## abd.set_current_model('ewens')
    ## print abd.theta, abd.I
    ## abd.set_current_model('etienne')
    ## print abd.theta, abd.I
    ## print '************************************************************\n'
    ## print 'LRT between Ewens and Etienne model (1 df): ',
    ## pv = abd.lrt('ewens', 'etienne')
    ## print pv, '\n'
    ## print 'set this model as current model:'
    ## abd.set_current_model('ewens' if pv > 0.05 else 'etienne')
    ## print abd
    ## print '************************************************************\n'
    ## t0 = time()
    ## gens = 500
    ## print 'Neutrality test p-value (under Ewens model):',
    ## print abd.test_neutrality(model='ewens', gens=gens)
    ## print 'should be arround: %s' % (0.9 if kind=='full' else 1.0)
    ## print '%s generations computed in %ss' % (gens, time()-t0)
    ## t0 = time()
    ## print 'Neutrality test p-value (under Etienne model):',
    ## print abd.test_neutrality(model='etienne', gens=gens)
    ## print 'Lognormality test p-value (under Etienne model):',
    ## abd.lognorm_optimal_params()
    ## print abd.test_neutrality(model='lognorm', gens=gens)
    ## ## print 'lognorm lnL', abd.lognorm_likelihood()
    ## ewens = abd.get_model('ewens')
    ## print 'ewens   lnL',ewens.lnL
    ## etienne = abd.get_model('etienne')
    ## print 'etienne lnL',etienne.lnL
    ## print 'should be arround: %s' % (0.1 if kind=='full' else 0.01)
    ## print '%s generations computed in %ss' % (gens, time()-t0)
    ## print '************************************************************\n'
    ## abd._kda = None
    ## #test_random_neutral(abd)
    ## if not test_failed:
    ##     print '\n\nAll test OK!\n'
    ## else:
    ##     print '\n\nFFFFailed in:', test_failed


if __name__ == "__main__":
    exit(main())
