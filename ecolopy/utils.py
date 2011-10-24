#!/usr/bin/python
"""
13 Jul 2011

some utils for abundance, essentially for computation of K(D,A)
"""

__author__  = "Francois-Jose Serra"
__email__   = "francois@barrabin.org"
__licence__ = "GPLv3"
__version__ = "0.12"

from gmpy2 import log, mul, mpfr, div, lngamma, gamma, exp
from sys import stdout
try:
    from matplotlib import pyplot
    import numpy as np
except ImportError:
    print 'matplotlib or numpt not found'

global STIRLINGS
STIRLINGS = {}


def generate_random_neutral_abundance (model_name, size, **kwargs):
    '''
    :argument model_name: model name (ewens, etienne, lognorm)
    :argument size: size of the community (J), if lognormal distribution size should be equal to number of species (S)
    :returns: random neutral distribution of abundance
    other args should be of kind theta, I, m
    
    **Example:**
    
    ::
    
      import ecolopy
      ecolopy.generate_random_neutral_abundance('ewens', 100, theta=12, I=12)
      
    '''
    # import inside because otherwise strange never-ending import...
    from ecolopy.ecological_model import EcologicalModel
    model = EcologicalModel (model_name, **kwargs)
    return model.rand_neutral(size)


def lpoch (z, m):
    '''
    returns log Pochhammer symbol taking advantage of:
    Pochhammer symbol (z)_m = (z)(z+1)....(z+m-1) = gamma(z+m)/gamma(z)
    '''
    return lngamma(z+m) - lngamma(z)


def shannon_entropy(abund, inds):
    '''
    computes Shannon entropy (H) for a given abundance table
    and number of individues.

    :argument abund: distribution of abundances as list
    :returns: shannon entropy
    '''
    return (sum ([-spe * log(spe) for spe in abund]) + inds * log (inds)) / inds


def table (out, spp=None):
    '''
    data to contingency table
    any kind of data

    :argument out: species abundance
    :returns: list of counts of individues by species
    '''
    setout = set (out)
    if spp == None:
        spp = len (setout)
    counts = dict (zip (setout, [mpfr(0.)]*spp))
    for ind in out:
        counts[ind] += mpfr(1)
    return [counts[x] for x in sorted (counts)]


def factorial_div (one, two):
    '''
    computes a!/b!
    '''
    if one < two:
        return div(1., reduce (mul, xrange(one, two)))
    elif one > two:
        return reduce (mul, xrange(two, one))
    else:
        return mpfr (1.)


def mul_polyn(polyn_a, polyn_b):
    '''
    computes the product of 2 polynomes, depending of the differences in length
    of the two polynomes, this function will call one of:
      * _mul_uneq_polyn: when length of polyn_a is >= length of polyn_b, will iterate over coeefficient.
      * _mul_simil_polyn: in case both polynomes have equal length, will iterate over indices.

    to test multiplication of pylnomes try equality of the two functions:
            _mul_uneq_polyn(polyn_a, polyn_b, len_a, len_b)
                                ==
            _mul_simil_polyn(polyn_a, polyn_b, len_a, len_b)

    **Example:**
    ::

      from ecolopy.utils import mul_polyn
      # (2 + 3^x + 5x^2)  * (x)
      mul_polyn([2,3,5], [0,1])
      # will return: [mpfr('0.0'), mpfr('2.0'), mpfr('3.0'), mpfr('5.0')]
      # that is: 0 + 2x + 3x^2 + 5x^3
      
    :argument polyn_a: list of indices of polynome
    :argument polyn_b: list of indices of polynome (e.g.: [1,3,5,0,2] for 2 + 3^x + 5x^2 + 0x^3 + 2x^4)
    :returns: a list representing multplication of the two polynomes
    '''
    if not polyn_a:
        return polyn_b
    if not polyn_b:
        return polyn_a
    # set default values
    len_a = len (polyn_a)
    len_b = len (polyn_b)
    diff = abs (len_a - len_b)
    if len_a >= len_b:
        polyn_b = polyn_b + [mpfr(0)] * (diff)
    else:
        _  = polyn_a + [mpfr(0.)] * (diff)
        polyn_a = polyn_b[:]
        polyn_b = _
        len_a = len (polyn_a)
        len_b = len (polyn_b)
    # switcher
    if len_a > len_b*2: # most
        return _mul_uneq_polyn(polyn_a, polyn_b, len_a, len_b)
    return _mul_simil_polyn(polyn_a, polyn_b, len_a, len_b)


def _mul_simil_polyn(polyn_a, polyn_b, len_a, len_b):
    '''
    fast polynomial multiplication when polynomes are nearly equal
    -> iterates over factors
    '''
    def mult_one (la, lb, stop, start=0):
        '''
        internal that computes row multiplication of 2 lists
        start and stops are here to skip multiplications by zero
        '''
        return [mul (la[i], lb[i]) for i in xrange (start, stop)]
    max_len = len_a + len_b
    diff = len_a - len_b
    new = []
    for i in xrange (1, len_a +1):
        new.append (sum (mult_one (polyn_a[:i], polyn_b[i-1::-1], i)))
    len_a2 = len_a * 2 - 1
    len_a1 = len_a + 1
    len_a3 = len_a - 1
    for i in xrange (len_a, max_len - 1):
        new.append (sum (mult_one (polyn_a[i-len_a3 : len_a1],
                              polyn_b[len_a    : i-len_a :-1], len_a2-i, diff)))
    return new    


def _mul_uneq_polyn(polyn_a, polyn_b, len_a, len_b):
    '''
    fast polynomial multiplication when 1 polynome >~ 2 times larger.
    -> iterates over coefficients
    '''
    new = [mpfr(0)] * (len_a + len_b - 1)
    for i in xrange (len_a):
        pai = polyn_a[i]
        for j in xrange (len_b):
            new [i + j] += pai * polyn_b[j]
    return new


def pre_get_stirlings(max_nm, needed, verbose=True):
    '''
    takes advantage of recurrence function:
      s(n, m) = s(n-1, m-1) - (n-1) * s(n-1, m)
    and as  s(whatever, 0) = 0 :
      s(n+1, 1) = -n * s(n, 1)
    keep only needed stirling numbers (necessary for large communities)

    :argument max_nm: max number of individues in a species
    :argument needed: list of individues count in species in our dataset, needed
                      in order to limit the number of stirling numbers kept in
                      memory.
    :argument True verbose: displays information about state.
    '''
    STIRLINGS [1, 1] = mpfr (1)
    for one in xrange (1, max_nm):
        if verbose and not one % 1000:
            stdout.write('\r    %s of %s, size: %s' % (one, max_nm,
                                                  STIRLINGS.__sizeof__()))
            stdout.flush()
        STIRLINGS [one+1, 1] = -one * STIRLINGS [one, 1]
        for two in xrange (2, one+1):
            STIRLINGS[one+1, two] = STIRLINGS[one, two-1] - one * STIRLINGS[one, two]
        STIRLINGS[one+1, one+1] = STIRLINGS[one, one]
        if one-1 in needed: continue
        for j in xrange (1, one):
            del (STIRLINGS[one-1,j])
    if verbose:
        stdout.write('\n')


def stirling (one, two):
    '''
    returns log unsingned stirling number, taking advantage of the fact that
    if x+y if odd signed stirling will be negative.
    takes also advantage of recurrence function:
    s(n,m) = s(n-1, m-1) - (n-1) * s(n-1,m)
    '''
    # return abs of stirling number
    if (one + two)%2:
        return -STIRLINGS [one, two]
    return STIRLINGS [one, two]


def fast_etienne_likelihood(abd, params, kda_x=None):
    '''
    same as Abundance inner function, but takes advantage of constant
    m value when varying only theta.

    :argument abd: Abundance object
    :argument params: list containing theta and m
    :argument kda_x: precimputed list of exp(kda + ind*immig)
    '''
    theta     = params[0]
    immig     = float (params[1]) / (1 - params[1]) * (abd.J - 1)
    log_immig = log (immig)
    theta_s   = theta + abd.S
    kda       = abd._kda
    if not kda_x:
        kda_x = [exp (kda [val] + val * log_immig) for val in xrange (abd.J - abd.S)]
    poch1 = exp (abd._factor + log (theta) * abd.S - \
                 lpoch (immig, abd.J) + \
                 log_immig * abd.S + lngamma(theta))
    gam_theta_s = gamma (theta_s)
    lik = mpfr(0.0)
    for val in xrange (abd.J - abd.S):
        lik += poch1 * kda_x[val] / gam_theta_s
        gam_theta_s *= theta_s + val
    return ((-log (lik)), kda_x)


def draw_contour_likelihood (abd, model=None, theta_range=None,
                             m_range=None, num_dots=100, local_optima=True,
                             write_lnl=False):
    """
    Draw contour plot of the log likelihood of a given abundance to fit Etienne model.   

    :argument abd: Abundance object
    :argument None model: model name, if None current model is used
    :argument None theta_range: minimum and maximum value of theta as list. If None, goes from 1 to number of species (S)
    :argument None m_range:  minimum and maximum value of m as list. If None, goes from 0 to 1
    :argument 100 num_dots: Number of dots to paint
    :argument True local_optima: display all local optima founds as white cross
    :argument False write_lnl: allow to write likelihood values manually by left click on the contour plot.
    """
    if not theta_range:
        theta_range = [1, int (abd.S)]
    if not m_range:
        m_range = [1e-16, 1-1e-9]
    if not model:
        model = abd.get_current_model_name()
    
    x = np.linspace(theta_range[0], theta_range[1], num_dots)
    y = np.linspace(m_range[0], m_range[1], num_dots)

    if model == 'etienne':
        lnl_fun = fast_etienne_likelihood
    elif model=='ewens':
        lnl_fun = lambda x: abd.ewens_likelihood(x[0])
    elif model == 'lognorm':
        lnl_fun = abd.lognorm_likelihood
    else:
        raise Exception ('Model not available.')

    abd.set_current_model(model)
    z = np.zeros ((num_dots, num_dots))
    if model=='etienne':
        kdas = {}
        for k, i in enumerate(x):
            print k, 'of', num_dots
            for l, j in enumerate (y):
                if not j in kdas:
                    lnl, kda_x = lnl_fun (abd, [i, j])
                    lnl = -lnl
                    kdas[j] = kda_x
                else:
                    lnl = -lnl_fun (abd, [i, j], kda_x=kdas[j])[0]
                z[l][k] = lnl
    else:
        for k, i in enumerate(x):
            print k, 'of', num_dots
            for l, j in enumerate (y):
                lnl = -lnl_fun ([i, j])
                z[l][k] = lnl

    if local_optima:
        loc_opt = []
        for i in xrange(1, len(z)-1):
            for j in xrange(1, len(z)-1):
                if z[i][j] > z[i][j-1] and \
                   z[i][j] > z[i][j+1] and \
                   z[i][j] > z[i-1][j] and \
                   z[i][j] > z[i+1][j]:
                    loc_opt.append ((i, j))
        for i, j in loc_opt:
            pyplot.plot(x[j], y[i], color='w', marker='x')

    # define number of color categories... perhaps not the best way
    perc10 = np.percentile(z, 10)
    perc90 = np.percentile(z, 90)
    levels = list (np.arange (perc90, perc10, -(perc90-perc10)/num_dots))[::-1] + \
             list (np.arange (perc90, z.max(), (z.max()-perc90)/num_dots))[1:] + [z.max()]
    pyplot.plot(abd.theta, abd.m, color='w', marker='*')
    pyplot.contourf (x, y, z, levels)
    pyplot.colorbar (format='%.3f')
    pyplot.vlines(abd.theta, 0, abd.m, color='w', linestyle='dashed')
    pyplot.hlines(abd.m, 0, abd.theta, color='w', linestyle='dashed')    
    if write_lnl:
        fig = pyplot.contour (x, y, z, levels)
        pyplot.clabel (fig, inline=1, fontsize=10, colors='black', manual=True)
    pyplot.title ("Log likelihood of abundance under Etienne model")
    pyplot.xlabel ("theta")
    pyplot.ylabel ("m")
    pyplot.axis (theta_range + m_range)
    pyplot.show()


def draw_shannon_distrib(neut_h, obs_h):
    '''
    draws distribution of shanon values for random neutral

    :argument neut_h: list of shannon entropies corresponding to simulation under neutral model
    :argument obs_h: shannon entropy of observed distribution of abundance
    
    '''
    neut_h = np.array ([float (x) for x in neut_h])
    obs_h = float (obs_h)
    pyplot.hist(neut_h, 40, color='green', histtype='bar', fill=True)
    pyplot.axvline(float(obs_h), 0, color='r', linestyle='dashed')
    pyplot.axvspan (float(obs_h) - neut_h.std(), float(obs_h) + neut_h.std(),
                    facecolor='orange', alpha=0.3)
    pyplot.xlabel('Shannon entropy (H)')
    pyplot.ylabel('Number of observations over %s simulations' % (len (neut_h)))
    pyplot.title("Histogram of entropies from %s simutations compared to \nobserved entropy (red), deviation computed from simulation" % (len (neut_h)))
    pyplot.show()


def mean (x):
    """
    :argument x: list of values
    :returns: mean of this list
    """
    return float (sum(x))/len (x)



def std (x):
    """
    :argument x: list of values
    :returns: standard deviation of this list
    """
    mean_x = mean (x)
    return (sum([(i - mean_x)**2 for i in x])/len (x))**.5
