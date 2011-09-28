#!/usr/bin/python
"""
19 Aug 2011


"""

__author__  = "Francois-Jose Serra"
__email__   = "francois@barrabin.org"
__licence__ = "GPLv3"
__version__ = "0.13"


from random import random
from utils  import table
from numpy  import exp
from scipy.stats import lognorm

def rand_neutral_etienne (inds, theta, immig):
    '''
    generates random distribution according to J, theta and I

    :arguments inds: number of individues in community (J)
    :arguments theta: corresponding to the model
    :arguments immig: immigration rate (I)

    :returns: distribution of abundance (list)
    
    '''
    theta = float (theta)
    immig = float (immig)
    mcnum   = [0] * int (inds)
    locnum  = [0] * int (inds)
    mcnum[0] = 1
    new = nxt = -1
    for ind in xrange (inds):
        if random () > immig / (ind + immig):
            locnum [ind] = locnum [int (random () * ind)]
        else:
            new += 1
            if random () <= theta / (theta + new):
                nxt += 1
                mcnum[new] = nxt + 1
            else:
                mcnum[new] = mcnum[int (random () * (new))]
            locnum[ind] = mcnum[new]
    return table (locnum, new + 1)


def rand_neutral_ewens (inds, theta):
    '''
    generates random distribution according to J and theta
    
    :arguments inds: number of individues in community (J)
    :arguments theta: corresponding to the model

    :returns: distribution of abundance (list)
    
    '''
    theta = float (theta)
    out = [0] * int (inds)
    out [0] = spp = 1
    for ind in xrange (inds):
        if random () < theta/(theta + ind):
            spp += 1
            out[ind] = spp
        else:
            out[ind] = out [int (random () * ind)]
    return table (out, spp)


def rand_lognormal (inds, sd, mu):
    '''
    generates random lognormal distribution
    
    :arguments inds: number of individues in community (J)
    :arguments sd: usually standard deviation of the distribution of abundaces
    :arguments mu: usually mean of the distribution of abundaces

    :returns: distribution of abundance (list)
    
    '''
    return [int (i+1) for i in lognorm.rvs (sd, scale=exp(mu),
                                            size=int (inds))]



