#!/usr/bin/python
"""
04 Jul 2011


"""

__author__  = "Francois-Jose Serra"
__email__   = "francois@barrabin.org"
__licence__ = "GPLv3"
__version__ = "0.1"

from scipy.stats    import chisqprob, lognorm
from scipy.optimize import fmin, fmin_slsqp, fmin_tnc, fmin_l_bfgs_b, golden
from gmpy2          import mpfr, log, exp, lngamma, gamma
from cPickle        import dump, load
from os.path        import isfile
from sys            import stdout
from numpy          import mean, std, ravel

from utils          import table, factorial_div, mul_polyn, shannon_entropy
from utils          import lpoch, lngamma, pre_get_stirlings, stirling
from utils          import del_stirling
from random_neutral import rand_neutral_etienne, rand_neutral_ewens
from random_neutral import rand_lognormal

class Abundance (object):
    '''
    Compute theta m I ect... for a given data sample
     * J_tot is the size of the metacommunity, if not defined = J * 3
    Values of theta, m, I etc will be stored in params dict, under key
    corresponding to model computed.
    '''

    def __init__ (self, data, j_tot=None):
        if type (data) != list:
            self.data_path = data
            data = self._parse_infile ()
        self.abund     = [mpfr(x) for x in sorted (data [:])]
        self.J         = mpfr(sum (data))
        self.S         = mpfr(len (data))
        self.j_tot     = j_tot if j_tot else self.J * 3
        self.shannon   = shannon_entropy (self.abund, self.J)
        self.params    = {}
        self._kda      = None


    def __repr__(self):
        """
        for print
        """
        return '''Abundance (object)
    Number of individues (J)  : %d
    Number of species (S)     : %d
    Shannon entropy (shannon) : %.4f
    Metacommunity size (j_tot): %d
    Models computed           : %s
    ''' % (self.J, self.S, self.shannon, self.j_tot,
           ', '.join (self.params.keys()))


    def lrt (self, model_1, model_2):
        '''
        likelihood-ratio-test between two neutral models
        returns p-value of rejection of alternative model
        '''
        return chisqprob(2 * (float (self.params[model_1]['lnL']) - \
                              float (self.params[model_2]['lnL'])), df=1)


    def ewens_likelihood (self, theta):
        '''
        get likelihood value of Ewens according to parthy/tetame
        returns likelihood
        '''
        factor = lngamma (self.J + 1)
        phi = table (self.abund)
        phi += [0] * int (max (self.abund) - len (phi))
        for spe in xrange (self.S):
            factor -= log (max (1, self.abund[spe]))
        for spe in xrange (max (self.abund)):
            factor -= lngamma (phi[spe] + 1)
        lnl = mpfr (lpoch (theta, self.J)) - log (theta) * self.S - factor
        self.params ['factor'] = factor
        return lnl

    def lognorm_likelihood (self, mu=None, sd=None):
        '''
        returns -log-likelihood of fitting to log-normal distribution
        '''
        data = [int (x) for x in self.abund]
        sd = std ([float (log (x)) for x in data]) if sd == None else sd
        if mu== None: mu = mean ([float (log (x)) for x in data])
        return -sum ([log (x) for x in lognorm.pdf (data, sd,
                                                    scale=float (exp (mu)))])

    def ewens_optimal_params (self):
        '''
        optimize theta using theta likelihood function, acording to Ewens
        model.
        '''
        self.params['ewens'] = tmp = {}
        theta_like   = lambda x: -self._ewens_theta_likelihood (x)
        tmp['theta'] = golden (theta_like, 
                                             brack=[.01/self.J, self.J])
        tmp['m']     = tmp['theta'] / self.J / mpfr(2)
        tmp['I']     = tmp['m'] * (self.J - 1) / (1 - tmp['m'])
        tmp['lnL']   = self.ewens_likelihood (tmp['theta'])


    def etienne_optimal_params (self, method='fmin'):
        '''
        optimize theta and I using etienne likelihood function
        using scipy package, values that are closest to the one proposed
        by tetame, are raised by fmin function:
          * fmin    : theta = 48.1838; I = 1088.19 ; lnL = -10091.166; t = 316s
        but lower likelihoods are found using other algorithms:
         -> bounds = [(1, self.S), (1e-50, 1 - 1e-50)]
          * slsqp   : theta = 143.43 ; I = 82.41   ; lnL = -10088.943; t = 227s
          * l_bfgs_b: theta = 140.11 ; I = 84.37   ; lnL = -10088.941; t = 343s
          * tnc     : theta = 52.62  ; I = 729.34  ; lnL = -10090.912; t = 297s
        '''
        # define bounds
        bounds = [(1, self.S), (1e-50, 1-1e-50)]
        # define starting values
        if 'ewens' in self.params:
            start = self.params ['ewens']['theta'], self.params['ewens']['m']
        else:
            start = self.S/2, 0.5
        # function minimization
        if   method == 'fmin':
            theta, mut = fmin (self.etienne_likelihood, start,
                               full_output=False, disp=0)
        elif method == 'slsqp':
            theta, mut  = fmin_slsqp (self.etienne_likelihood, start,
                                      bounds=bounds,iprint=9)
        elif method == 'l_bfgs_b':
            theta, mut = fmin_l_bfgs_b (self.etienne_likelihood, start,
                                       bounds=bounds, approx_grad=True)[0]
        elif method == 'tnc':
            theta, mut = fmin_tnc (self.etienne_likelihood, start,
                                  bounds=bounds, approx_grad=True)[0]
        self.params['etienne']['theta'] = theta
        self.params['etienne']['m']     = mut
        self.params['etienne']['I']     = mut * (self.J - 1) / (1 - mut)
        self.params['etienne']['lnL']   = self.etienne_likelihood ([theta, mut])


    def _ewens_theta_likelihood (self, theta):
        '''
        returns the likelihood of theta for a given dataset
        '''
        if theta < 0:
            return mpfr ('-inf')
        return self.S * log (theta) + lngamma (theta) - lngamma (theta + self.J)


    def test_neutrality (self, model='ewens', gens=100):
        '''
        test for neutrality comparing shanon entropy
        returns p_value
        '''
        pval = 0
        if model == 'lognorm':
            logabund = [float (log (x)) for x in self.abund]
            theta = mean (logabund)
            immig = std  (logabund)
            print theta, immig
        else:
            theta = self.params[model]['theta']
            immig = self.params[model]['I']
        for _ in xrange (gens):
            pval += shannon_entropy (self.rand_neutral (theta, immig,
                                                        model=model),
                                     self.J) < self.shannon
        return float (pval)/gens
    

    def rand_neutral (self, theta, immig, model='etienne'):
        '''
        Depending on model, generates random neutral abundance
        with theta and I
        '''
        if model == 'ewens':
            return rand_neutral_ewens (self.J, theta)
        elif model == 'etienne':
            return rand_neutral_etienne (self.J, theta, immig)
        elif model == 'lognorm':
            return rand_lognormal (self.S, immig, theta)


    def _parse_infile (self):
        '''
        parse infile and return list of abundances
        '''
        abundances = []
        for line in open (self.data_path):
            abundances.append (int (line.strip()))
        return abundances


    def etienne_likelihood(self, params, verbose=True):
        '''
        logikelihood function where
          params[0] = theta
          params[1] = I
          K[A]      = log (K(D,A))
          K[A]      = K(D,A) in Etienne paper
        cf Etienne, 2005
        x = [48.1838, 1088.19]
        returns log likelihood of given theta and I
        '''
        if not 'factor' in self.params: # define it
            self.ewens_optimal_params()
        if not 'etienne' in self.params:
            if verbose:
                print "\nGetting K(D,A) according to Etienne 2005 formula:"
            self._get_kda (verbose=verbose)
        kda       = self._kda
        theta     = params[0]
        immig     = float(params[1])/(1 - params[1]) * (self.J - 1)
        log_immig = log (immig)
        theta_s   = theta + self.S
        poch1 = exp (self.params['factor'] + log (theta) * self.S - \
                     lpoch (immig, self.J) + \
                     log_immig * self.S + lngamma(theta))
        gam_theta_s = gamma (theta_s)
        lik = mpfr(0.0)
        for abd in xrange (self.J - self.S):
            lik += poch1 * exp (kda [abd] + abd * log_immig) / gam_theta_s
            gam_theta_s *= theta_s + abd
        return -log (lik)


    def dump_params (self, outfile):
        '''
        save params and kda with pickle
        '''
        if isfile (outfile):
            self.load_params (outfile)
        self.params['KDA'] = self._kda[:]
        dump (self.params, open (outfile, 'w'))
        del (self.params['KDA'])


    def load_params (self, infile):
        '''
        load params and kda with pickle from infile
        '''
        self.params.update (load (open (infile)))
        if 'KDA' in self.params:
            self._kda = self.params['KDA']
            del (self.params['KDA'])


    def _get_kda (self, verbose=True):
        '''
        compute kda according to etienne formula
        '''
        specabund = [sorted (list (set (self.abund))), table (self.abund)]
        sdiff     = len (specabund [1])
        polyn     = []
        # compute all stirling numbers taking advantage of recurrence function
        needed = {0: True}
        for i in xrange (sdiff):
            for k in xrange (1, specabund[0][i] + 1):
                needed [int (specabund[0][i])] = True
        if verbose:
            stdout.write('  Getting some stirling numbers...\n')
        pre_get_stirlings (max (specabund[0]), needed)
        for i in xrange (sdiff):
            if verbose:
                stdout.write ("\r  Computing species %s out of %s" % (i+1,
                                                                      sdiff))
                stdout.flush ()
            polyn1 = []
            for k in xrange (1, specabund[0][i] + 1):
                coeff = stirling (int (specabund[0][i]), k) * \
                        factorial_div (k, specabund[0][i])
                polyn1.append (coeff)
            if not polyn1:
                polyn1.append(mpfr(1.))
            # polyn1 exponential the number of individues for current species
            polyn2 = polyn1[:]
            for _ in xrange (1, specabund[1][i]):
                polyn1 = mul_polyn (polyn1, polyn2)
            # polyn = polyn * polyn1
            polyn = mul_polyn (polyn, polyn1)
        kda = []
        stdout.write ('\n')
        for i in polyn:
            kda.append (log (i))
        self.params.setdefault ('etienne', {})
        self._kda = kda
        #del_stirling()


