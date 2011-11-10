#!/usr/bin/python
"""
04 Jul 2011


"""

__author__  = "Francois-Jose Serra"
__email__   = "francois@barrabin.org"
__licence__ = "GPLv3"
__version__ = "0.13"

from scipy.stats      import chisqprob, lognorm
from scipy.optimize   import fmin, fmin_slsqp, fmin_tnc, fmin_l_bfgs_b, golden
from gmpy2            import mpfr, log, exp, lngamma, gamma
from cPickle          import dump, load
from os.path          import isfile
from sys              import stdout, stderr
                      
from utils    import table, factorial_div, mul_polyn, shannon_entropy
from utils    import lpoch, pre_get_stirlings, stirling## , mean , std
from ecolopy.ecological_model import EcologicalModel

class Abundance (object):
    '''
    Main abundance class.

    :argument data: containing species count, can be given in one of those format:
    
       * python list, each element being a species count
       * text file containing one species count per line
       * pickle file containing abundance object

    :argument None j_tot: is the size of the metacommunity, if not defined j_tot = J * 3

    :returns: an abundance object

    **Examples:**
    ::

      # python list
      abd = Abundance([1, 1, 2, 3, 4, 8, 12])
      # text file
      abd = Abundance("bci.txt", j_tot=256987)
      # finally a saved Abundance object
      abd = Abundance("abundance_bci_first_try.pik")
    
    '''

    def __init__ (self, data, j_tot=None):
        """
        creation of Abundance object
        """
        if type (data) != list:
            self.data_path = data
            data = self._parse_infile ()
        if not data is None:
            self.abund     = [mpfr(x) for x in sorted (data [:])]
            self._kda      = None
        self.J         = mpfr(sum (self.abund))
        self.S         = mpfr(len (self.abund))
        self.j_tot     = j_tot if j_tot else self.J * 3
        self.shannon   = shannon_entropy (self.abund, self.J)
        self.theta     = None
        self.m         = None
        self.I         = None
        self.lnL       = None
        self._factor   = None
        self.__models  = {}
        self.__current_model = None


    def __str__(self):
        """
        for print
        """
        return '''Abundance (object)
    Number of individuals (J) : %d
    Number of species (S)     : %d
    Shannon entropy (shannon) : %.4f
    Metacommunity size (j_tot): %d
    Models computed           : %s
    Current model             : %s
    Theta                     : %s
    I                         : %s
    m                         : %s
    lnL                       : %s
    ''' % (self.J, self.S, self.shannon, self.j_tot,
           ', '.join (self.__models.keys()),
           self.__current_model.name if self.__current_model else None,
           self.theta, self.I, self.m, self.lnL)


    def lrt (self, model_1, model_2, df=1):
        '''
        likelihood-ratio-test between two neutral models

        :argument model_1: string representing simplest model, usually Ewens
        :argument model_2: string representing most complex model, usually Etienne
        :argument 1 df: number of degrees of freedom (1 when comparing Etienne and Ewens)
        
        :returns: p-value of rejection of alternative model
        
        eg: usually ewens, etienne. And if pval < 0.05 than use etienne
        '''
        return chisqprob(2 * (float (self.get_model(model_1).lnL) - \
                              float (self.get_model(model_2).lnL)), df=df)


    def get_current_model_name (self):
        """

        :returns: current model name
        
        """
        return self.__current_model.name

    def iter_models (self):
        """
        iterate over pre computed models

        :returns: model object
        
        """
        for model in self.__models.values():
            yield model


    def set_model (self, model):
        """
        add one model computed externally to the computed models of current Abundance object

        :argument model: model object
        
        """
        self.__models[model.name] = model


    def get_model (self, name):
        """
        :argument name: name of a computed model
        :returns: a EcologicalModel object corresponding to one of the computed models
        """
        if name in self.__models:
            return self.__models[name]
        else:
            return None

    def ewens_likelihood (self, theta):
        '''
        get likelihood value of Ewens according to parthy/tetame
        
        :argument theta: value of theta
        :returns: likelihood
        '''
        factor = lngamma (self.J + 1)
        phi = table (self.abund)
        phi += [0] * int (max (self.abund) - len (phi))
        for spe in xrange (self.S):
            factor -= log (max (1, self.abund[spe]))
        for spe in xrange (max (self.abund)):
            factor -= lngamma (phi[spe] + 1)
        lnl = mpfr (lpoch (theta, self.J)) - log (theta) * self.S - factor
        self._factor = factor
        return lnl


    ### NO LONGER USED
    ## def lognorm_likelihood (self, params=None):
    ##     '''
    ##     returns -log-likelihood of fitting to log-normal distribution
    ## 
    ##     :argument None params: a list of 2 parameters:
    ##     
    ##       * mu = parmas[0]
    ##       * sd = parmas[1]
    ##     :returns: log likelihood given the parameters
    ##     
    ##     '''
    ##     data = [int (x) for x in self.abund]
    ##     if params == None:
    ##         log_data = [float (log (x)) for x in data]
    ##         mu, sd = mean (log_data), std (log_data)
    ##     else:
    ##         mu, sd = params
    ##     return -sum ([log (x) for x in lognorm.pdf (data, sd,
    ##                                                 scale=float (exp (mu)))])


    def etienne_likelihood(self, params):
        '''
        logikelihood function
        
        :argument params: a list of 2 parameters:
        
          * theta = parmas[0]
          * I     = parmas[1]
        :returns: log likelihood of given theta and I
        
        '''
        if not self._factor: # define it
            self.ewens_optimal_params()
        kda       = self._kda
        theta     = params[0]
        immig     = float (params[1]) / (1 - params[1]) * (self.J - 1)
        log_immig = log (immig)
        theta_s   = theta + self.S
        poch1 = exp (self._factor + log (theta) * self.S - \
                     lpoch (immig, self.J) + \
                     log_immig * self.S + lngamma(theta))
        gam_theta_s = gamma (theta_s)
        lik = mpfr(0.0)
        for abd in xrange (self.J - self.S):
            lik += poch1 * exp (kda [abd] + abd * log_immig) / gam_theta_s
            gam_theta_s *= theta_s + abd
        return -log (lik)


    def ewens_optimal_params (self):
        '''
        Main function to optimize theta using theta likelihood function, acording to Ewens
        model.
        '''
        theta_like   = lambda x: -self._ewens_theta_likelihood (x)
        theta = golden (theta_like, 
                                        brack=[.01/self.J, self.J])
        m     = theta / self.J / mpfr(2)
        I     = m * (self.J - 1) / (1 - m)
        lnL   = self.ewens_likelihood (theta)
        self.__models['ewens'] = EcologicalModel ('ewens', theta=theta,
                                                  m=m, I=I, lnL=lnL)


    def lognorm_optimal_params (self):
        '''
        Main function to get optimal params for lognormal distribution according to abundance
        '''
        ## USED BEFORE WITH LOGNORM_LIKELIHOOD
        ## data = self.abund
        ## start = (mean ([float (log (x)) for x in data]),
        ##          std ([float (log (x)) for x in data]))
        ## mu, sd = fmin (self.lognorm_likelihood, start,
        ##                full_output=False, disp=0)
        ## lnL   = self.lognorm_likelihood ((mu, sd))
        shape, loc, scale = lognorm.fit([float (x) for x in self.abund])
        self.__models['lognorm'] = EcologicalModel ('lognorm', theta=shape,
                                                    I=loc, m=scale, lnL=None)


    def etienne_optimal_params (self, method='fmin', verbose=True):
        '''
        Main function to optimize theta and I using etienne likelihood function
        using scipy package, values that are closest to the one proposed
        by tetame, are raised by fmin function.

        :argument fmin method: optimization strategy, can be one of fmin, slsqp, l_bfgs_b or tnc (see scipy.optimize documentation)
        
        '''
        # define bounds
        bounds = [(1, self.S), (1e-49, 1-1e-49)]
        all_ok  = True
        # define starting values
        if 'ewens' in self.__models:
            m = self.get_model ('ewens')
            start = m.theta, m.m
        else:
            start = self.S/2, 0.5
        # conpute kda
        if not 'etienne' in self.__models:
            if verbose:
                print "\nGetting K(D,A) according to Etienne 2005 formula:"
            self._get_kda (verbose=verbose)
        # function minimization
        if   method == 'fmin':
            theta, mut = fmin (self.etienne_likelihood, start,
                               full_output=False, disp=0)
        elif method == 'slsqp':
            a, _, _, err, _ = fmin_slsqp (self.etienne_likelihood, start,
                                          iter=1000, iprint=0,
                                          bounds=bounds, full_output=True)
            theta, mut = a
            if err != 0:
                all_ok = False
        elif method == 'l_bfgs_b':
            a, _, err = fmin_l_bfgs_b (self.etienne_likelihood, start,
                                       maxfun=1000, bounds=bounds,
                                       iprint=-1, approx_grad=True)
            theta, mut = a
            if err['warnflag'] != 0:
                all_ok = False
        elif method == 'tnc':
            a, _, err = fmin_tnc (self.etienne_likelihood, start, maxfun=1000,
                                  messages=0, bounds=bounds,
                                  approx_grad=True)
            theta, mut = a
            if err != 1:
                all_ok = False
        I     = mut * (self.J - 1) / (1 - mut)
        lnL   = self.etienne_likelihood ([theta, mut])
        self.__models['etienne'] = EcologicalModel ('etienne', theta=theta,
                                                    m=mut, I=I, lnL=lnL)
        if not all_ok:
            raise Exception ('Optimization failed')


    def set_current_model (self, name):
        """
        set on model as default/current model.

        :argument name: model name of precomputed model
        
        """
        self.__current_model = self.get_model(name)
        for key in ['theta', 'I', 'm', 'lnL']:
            try:
                self.__dict__[key] = self.__current_model.__dict__[key]
            except KeyError:
                self.__dict__[key] = None
                

    def _ewens_theta_likelihood (self, theta):
        '''
        returns the likelihood of theta for a given dataset
        
        '''
        if theta < 0:
            return mpfr ('-inf')
        return self.S * log (theta) + lngamma (theta) - lngamma (theta + self.J)


    def test_neutrality (self, model='ewens', gens=100, give_h=False, fix_s=False, tries=1000):
        '''
        test for neutrality comparing shanon entropy
        if (Hobs > Hrand_neut) then eveness of observed data is
        higher then neutral

        :argument ewens model: model name otherwise, current model is used
        :argument 100 gens: number of random neutral distributions to generate
        :argument False give_h: also return list of shannon entropies
        :argument False fix_s: decide wether to fix or not the number of species for the generation of random neutral communities
        :argument False tries: in case S is fixed, determines the number of tries in order to obtain the exact same number of species as original comunity. In case The number of tries is exceeded, an ERROR message is displayed, and p-value returned is 1.0.
        :returns: p_value anf if give_h also returns shannon entropy of
        all random neutral abundances generated
        '''
        fast_shannon = lambda abund: sum ([-spe * log(spe) for spe in abund])
        pval = 0
        inds = self.S if model is 'lognorm' else self.J
        model = self.get_model (model)
        if not model:
            return None
        neut_h = []
        for _ in xrange (gens):
            if fix_s:
                for _ in xrange (tries):
                    tmp = model.rand_neutral (inds)
                    if len(tmp) == self.S:
                        break
                else:
                    stderr.write('ERROR: Unable to obtain S by simulation')
                    return 1.0
            else:
                tmp = model.rand_neutral (inds)
            l_tmp = sum (tmp)
            neut_h.append ((fast_shannon (tmp) + l_tmp*log(l_tmp))/l_tmp)
            pval += neut_h[-1] < self.shannon
        if give_h:
            return float (pval)/gens, neut_h
        return float (pval)/gens
    

    def _parse_infile (self):
        '''
        parse infile and return list of abundances
        '''
        abundances = []
        lines = open (self.data_path).readlines()
        if lines[0].strip() == '(dp1':
            self.load_abundance (self.data_path)
            return None
        for line in lines:
            abundances.append (int (line.strip()))
        return abundances



    def dump_abundance (self, outfile, force=False):
        '''
        save params and kda with pickle
        force option is for writing pickle with no consideration
        if existing

        :argument outfile: path of the outfile
        :argument False force: overwrite existing outfile
        
        '''
        if isfile (outfile) and not force:
            self.load_abundance (outfile)
        try:
            self.__models['KDA']   = self._kda[:]
        except TypeError:
            self.__models['KDA'] = None
        self.__models['ABUND'] = self.abund[:]
        dump (self.__models, open (outfile, 'w'))
        del (self.__models['KDA'])
        del (self.__models['ABUND'])

    def load_abundance (self, infile):
        '''
        load params and kda with pickle from infile
        WARNING: do not overright values of params.

        :argument infile: path of the outfile

        '''
        old_params = load (open (infile))
        old_params.update(self.__models)
        self.__models = old_params
        self.abund = self.__models['ABUND']
        del (self.__models['ABUND'])
        if 'KDA' in self.__models:
            self._kda  = self.__models['KDA']
            del (self.__models['KDA'])
        else:
            self._kda = None


    def generate_random_neutral_distribution (self, model=None, J=None):
        """
        return distribution of abundance, according to a given model
        and a given community size.
        if none of them are given, values of current Abundance are used.

        :argument None model: model name (default current model)
        :argument None J: size of wanted community (default size of the community of current abundance)
        :returns: random neutral distribution of abundances
        
        """
        if not model:
            try:
                model = self.__current_model
            except:
                return None
        else:
            model = self.get_model(model)
        if model is None:
            raise Exception ('Need first to optimize this model,\n       ' + \
                             '    Unable to generate random distribution.')
        if not J:
            J = self.S if model.name is 'lognorm' else self.J
        return model.rand_neutral (J)

    
    def _get_kda (self, verbose=True):
        '''
        compute K(D,A) according to etienne formula
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
        pre_get_stirlings (max (specabund[0]), needed, verbose=verbose)
        for i in xrange (sdiff):
            if verbose:
                stdout.write ("\r  Computing K(D,A) at species %s out of %s" \
                              % (i+1, sdiff))
                stdout.flush ()
            polyn1 = []
            for k in xrange (1, specabund[0][i] + 1):
                coeff = stirling (specabund[0][i], k) * \
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
        if verbose:
            stdout.write ('\n')
        for i in polyn:
            kda.append (log (i))
        self._kda = kda



