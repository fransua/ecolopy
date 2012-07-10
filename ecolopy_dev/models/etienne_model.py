#!/usr/bin/python
"""
21 Sep 2011


"""

__author__  = "Francois-Jose Serra"
__email__   = "francois@barrabin.org"
__licence__ = "GPLv3"
__version__ = "0.13"


from ecolopy_dev.models.untb_model import UNTBModel
from ecolopy_dev.utils             import table_mpfr, table, factorial_div
from ecolopy_dev.utils             import power_polyn, lpoch
from ecolopy_dev.utils             import pre_get_stirlings, stirling, mul_polyn
from warnings                      import warn

try:
    from gmpy2                         import log, lngamma, exp, gamma, mpfr
except ImportError:
    warn("WARNING: GMPY2 library not found, using numpy")
    from numpy                         import log, exp, float128 as mpfr
    from scipy.special                 import gamma, gammaln as lngamma
    
from scipy.optimize                import fmin, fmin_slsqp, fmin_tnc, fmin_l_bfgs_b
from sys                           import stdout
from random                        import random


class EtienneModel(UNTBModel):
    """
    Class representing Ecological models

    :argument name: name of the class, can be either ewens, etienne or lognorm
    :returns: EcologicalModel object
    """
    def __init__(self, community, **kwargs):
        if not 'kda' in kwargs:
            self._kda = None
        super(EtienneModel, self).__init__(community, **kwargs)
        self.optimize(**kwargs)


    def random_community (self, inds=None, theta=None, immig=None):
        '''
        generates random distribution according to J, theta and I

        :argument inds: number of individuals in community (J)
        :argument theta: corresponding to the model
        :argument immig: immigration rate (I)

        :returns: distribution of abundance (list)

        '''
        theta    = float (theta) if theta else self.theta
        inds     = inds or self.community.J
        immig    = float (immig) if immig else self.I
        mcnum    = [0] * int (inds)
        locnum   = [0] * int (inds)
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


    def likelihood(self, params):
        '''
        log-likelihood function
        
        :argument params: a list of 2 parameters:
        
          * theta = params[0]
          * m     = params[1]
        :returns: log likelihood of given theta and I
        
        '''
        kda       = self._kda
        theta     = params[0]
        immig     = float (params[1]) / (1 - params[1]) * (self.community.J - 1)
        log_immig = log (immig)
        theta_s   = theta + self.community.S
        poch1 = exp (self._factor + log (theta) * self.community.S - \
                     lpoch (immig, self.community.J) + \
                     log_immig * self.community.S + lngamma(theta))
        gam_theta_s = gamma (theta_s)
        lik = mpfr(0.0)
        for abd in xrange (self.community.J - self.community.S):
            lik += poch1 * exp (kda [abd] + abd * log_immig) / gam_theta_s
            gam_theta_s *= theta_s + abd
        return -log (lik)

        
    def optimize (self, method='fmin', start=None, verbose=True):
        '''
        Main function to optimize theta and I using etienne likelihood function
        using Scipy package, values that are closest to the one proposed
        by Tetame, are raised by fmin function.

        :argument fmin method: optimization strategy, can be one of fmin, slsqp, l_bfgs_b or tnc (see scipy.optimize documentation)
        :argument (community.S,0.5) start: tupple for startin values of theta and m
        :argument True verbose: displays running status
        
        '''
        # define bounds
        bounds = [(0, len (self.community.abund)), (1e-49, 1-1e-49)]
        all_ok  = True
        err = ''
        # define starting values
        start = start or self.community.S/2, 0.5
        # compute kda
        if not self._kda:
            if verbose:
                print "\nGetting K(D,A) according to Etienne 2005 formula:"
            self._get_kda (verbose=verbose)
        # function minimization
        if   method == 'fmin':
            theta, mut = fmin (self.likelihood, start,
                               full_output=False, disp=0)
        elif method == 'slsqp':
            a, _, _, err, _ = fmin_slsqp (self.likelihood, start,
                                          iter=1000, iprint=0,
                                          bounds=bounds, full_output=True)
            theta, mut = a
            if err != 0:
                all_ok = False
        elif method == 'l_bfgs_b':
            a, _, err = fmin_l_bfgs_b (self.likelihood, start,
                                       maxfun=1000, bounds=bounds,
                                       iprint=-1, approx_grad=True)
            theta, mut = a
            if err['warnflag'] != 0:
                all_ok = False
        elif method == 'tnc':
            a, _, err = fmin_tnc (self.likelihood, start, maxfun=1000,
                                  messages=0, bounds=bounds,
                                  approx_grad=True)
            theta, mut = a
            if err != 1:
                all_ok = False
        self._parameters['I']     = mut * (self.community.J - 1) / (1 - mut)
        self._lnL                 = self.likelihood ([theta, mut])
        self._parameters['m']     = mut
        self._parameters['theta'] = theta
        if not all_ok:
            raise Exception ('Optimization failed\n', err)


    def _get_kda (self, verbose=True):
        '''
        compute K(D,A) according to etienne formula
        '''
        specabund = [sorted (list (set (self.community.abund))),
                     table_mpfr (self.community.abund)]
        sdiff     = len (specabund [1])
        polyn     = []
        # compute all stirling numbers taking advantage of recurrence function
        needed = {0: True}
        for i in xrange (sdiff):
            for k in xrange (1, specabund[0][i] + 1):
                needed [int (specabund[0][i])] = True
        if verbose:
            stdout.write('  Getting some stirling numbers...\n')
        pre_get_stirlings (max (specabund[0]), needed, verbose=verbose)
        for i in xrange (sdiff):
            if verbose:
                stdout.write ("\r  Computing K(D,A) at species %s out of %s" \
                              % (i+1, sdiff))
                stdout.flush ()
            sai = specabund[0][i]
            polyn1 = [stirling(sai, k)*factorial_div(k, sai) for k in xrange (1, sai+1)]
            polyn1 = power_polyn (polyn1, specabund[1][i]) if polyn1 else [mpfr(1.)]
            polyn = mul_polyn (polyn, polyn1)
        self._kda = [log (i) for i in polyn]
        if verbose:
            stdout.write ('\n')

