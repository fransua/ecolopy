#!/usr/bin/python
"""
02 Jan 2012


"""

__author__  = "Francois-Jose Serra"
__email__   = "francois@barrabin.org"
__licence__ = "GPLv3"
__version__ = "0.0"


from random import random, lognormvariate
import numpy as np
from numpy.random import geometric
from sys import argv
from lockfile import FileLock
from cPickle import load, dump

from ecolopy_dev.community import Community
from ecolopy_dev.random_neutral import rand_neutral_ewens, rand_neutral_etienne

def get_random_S():
    return int (10+random()*990)

def generate_lognorm_distribution(S, mu=None, sd=None):
    if not mu:
        mu = 2+random()*2
    if not sd:
        sd = 1+random()*2
    return [lognormvariate(mu, sd) for _ in xrange (S)]

    
def generate_geometric_distribution(S,pr=None):
    if not pr:
        pr = random()
    return ([geometric(pr) for _ in xrange(S)])
    

def generate_ewens_distribution(S, inds):
    of = False
    for _ in xrange (1000):
        theta = random()*S
        if inds == 0: continue
        abd = []
        for j in xrange(100):
            try:
                abd = rand_neutral_ewens(inds, theta)
            except IndexError, e:
                print j, theta, inds
                raise e
            if len (abd) == S:
                of = True
                break
        if of:
            break
    else:
        raise Exception('failed')
    return abd


def generate_etienne_distribution (S, inds):
    of = False
    for _ in xrange(1000):
        theta = random()*S
        m = random()
        immig = m*(inds-1)/(1-m)
        for j in xrange(100):
            try:
                abd = rand_neutral_etienne(inds, theta, immig)
            except IndexError, e:
                print j, theta, immig, m, inds
                raise e
            if len (abd) == S:
                of = True
                break
        if of:
            break
    else:
        raise Exception('failed')
    return abd



def generate_random_distribution(model, S, inds):
    if model == 'etienne':
        distr = generate_etienne_distribution(S, inds)
    elif model == 'ewens':
        distr = generate_ewens_distribution(S, inds)
    elif model == 'lognorm':
        distr = generate_lognorm_distribution(S)
    elif model == 'geometric':
        distr = generate_geometric_distribution(S)
    return [int (float(i)+0.5) for i in distr if float(i)>=0.5]


def main():
    """
    main function
    """

    sizes = [10]#,30,50]#,100,150]#,200,400,1000]
    people  = [500, 1000]#100000]#, 1000000]
    models = ['geometric','lognorm','ewens','etienne']
    results = {}
    for typ in ['fix','lib']:
        results[typ] = {}
        for m1 in models:
            results[typ][m1]  = {}
            results[typ][m1]['all']  = []
            for S in sizes:
                results[typ][m1][S] = {}
                for inds in people + ['all']:
                    if inds == 1000000 and S < 50:
                        continue
                    results[typ][m1][S][inds] = []
    for m in models:
        for S in sizes:
            for inds in people:
                if inds == 1000000 and S < 50:
                    continue
                print m
                distr = generate_random_distribution(m, S, inds)
                print 'S: {0}, J: {1}'.format (len (distr), sum(distr))
                abd = Community(distr)
                abd.fit_model('ewens')
                abd.fit_model('etienne', verbose=False)
                lrt = abd.lrt('ewens', 'etienne')
                best = 'ewens' if lrt > 0.05 else 'etienne'
                pv_fs = abd.test_neutrality(gens=1000, model=best, n_cpus=8)
                pv_ls = abd.test_neutrality(gens=1000, model=best, fix_s=True, n_cpus=8)
                pv_fl = abd.test_neutrality(gens=1000, model=best,
                                            method='loglike', n_cpus=8)
                pv_ll = abd.test_neutrality(gens=1000, model=best, fix_s=True,
                                           method='loglike', n_cpus=8)
                results['fix'][m][S][inds].append(pv_fs)
                results['lib'][m][S][inds].append(pv_ls)
                results['fix'][m][S]['all'].append(pv_fs)
                results['lib'][m][S]['all'].append(pv_ls)
                results['fix'][m]['all'].append(pv_fs)
                results['lib'][m]['all'].append(pv_ls)
                

    log = '/tmp/abd_neut_test_simulations_01.pik'
    lock = FileLock(log)
    print log
    with lock:
        try:
            old_results = load(open(log))
        except:
            old_results = {}
        for typ in ['fix','lib']:
            old_results.setdefault(typ, {})
            for m1 in models:
                old_results[typ].setdefault(m1, {})
                for S in sizes:
                    old_results[typ][m1].setdefault(S, {})
                    for inds in people + ['all']:
                        if inds == 1000000 and S<50:
                            continue
                        old_results[typ][m1][S].setdefault(inds, [])
                        old_results[typ][m1][S][inds].append(results[typ][m1][S][inds])
                old_results[typ][m1].setdefault('all', [])
                old_results[typ][m1]['all'].append(results[typ][m1]['all'])
        pik = open(log, 'w')
        dump(old_results, pik)
        pik.close()

exit(main())



