#!/usr/bin/python
"""
04 Jul 2011


"""

__author__  = "Francois-Jose Serra"
__email__   = "francois@barrabin.org"
__licence__ = "GPLv3"
__version__ = "0.0"


from scipy import optimize
from numpy import log
from scipy.special import gammaln
from random import random
from scipy.stats import chisqprob

def table (out, spp=None):
    '''
    data to contingency table
    any kind of data
    '''
    if spp == None:
        spp = max (out)
    counts = dict (zip (set (out), [0]*spp))
    for ind in out:
        counts[ind] += 1
    return counts.values()

lpochham = lambda x, n: gammaln(x+n)-gammaln(x)

class Abundance (object):
    '''
    Compute theta m I ect... for a given data sample
    '''
    def __init__ (self, data, J_tot=None):
        if type (data) != list:
            self.data_path = data
            data = self.parse_infile ()
        self.abund = data [:]
        self.J     = sum (data)
        self.S     = len (data)
        self.theta = self._ewens_optimal_theta ()
        self.m     = self.theta / self.J / 2
        self.I     = self.m * float (self.J - 1) / (1 - self.m)
        self.J_tot = J_tot if J_tot else self.J * 3
        self.H     = self._shannon_entropy ()

    def _lrt_ewens_etienne(self):
        '''
        likelihood-ratio-test between two neutral models
        returns p-value of rejection of alternative model
        '''
        chisqprob(2 * (self.ewens_likelihood() - 0), df=1)

    def ewens_likelihood (self):
        '''
        get likelihood value of Ewens according to parthy/tetame
        returns likelihood and factor
        '''
        factor = gammaln (self.J+1)
        phi = table (self.abund)
        phi += [0] * (max (self.abund)-len (phi))
        for spe in xrange (self.S):
            factor -= log (max (1, self.abund[spe]))
        for spe in xrange (max (self.abund)):
            factor -= gammaln (phi[spe] + 1)
        return (lpochham (self.theta,
                          self.J) - log (self.theta) * self.S - factor,
                factor)

    def _ewens_optimal_theta (self):
        '''
        optimize theta using theta likelihood function
        '''
        theta_like = lambda x: -self._ewens_theta_likelihood (x)
        return optimize.golden(theta_like, brack=[.01/self.J, self.J])

    def _ewens_theta_likelihood (self, theta):
        '''
        returns the likelihood of theta for a given dataset
        '''
        theta = float (theta)
        return self.S * log (theta) + gammaln (theta) - gammaln (theta + self.J)

    def _shannon_entropy (self):
        '''
        computes shannon entropy
        '''
        H = 0
        for spe in self.abund:
            H -= spe * log (spe)
        return (H + self.J * log(self.J)) / self.J

    def test_neutrality (self, gens=100):
        '''
        test for neutrality comparing shanon entropy
        returns p_value
        '''
        def local_shannon (data, J):
            '''
            same as self.shannon entropy
            '''
            H = 0
            for spe in data:
                H -= spe * log (spe)
            return (H + J * log(J)) / J
        pval = 0
        for _ in xrange (gens):
            pval += local_shannon (self.rand_neutral(), self.J) < self.H
        return float (pval)/gens
    
    def rand_neutral (self):
        '''
        J should be the number of individuals in the dataset, and theta,
        the optimal theta of the dataset accod
        J = 376
        theta = 9.99
        '''
        theta = float (self.theta)
        out = [0] * self.J
        out [0] = spp = 1
        for ind in xrange (1, self.J):
            if random () < theta/(theta + ind):
                spp += 1
                out [ind] = spp
            else:
                out [ind] = out [int (random () * ind)]
        return table (out, spp)
        
    def parse_infile (self):
        '''
        parse infile and return list of abundances
        '''
        abundances = []
        for line in open (self.data_path):
            abundances.append (int (line.strip()))
        return abundances


    def etienne_likelihood(self, x):
        '''
        logikelihood function where
          x[0]=theta
          x[1]=I,K[A]=log(K(D,A))
          
        K[A]=K(D,A)in Etienne paper
        
        cf Etienne, 2005
        '''
        divisor = 0.0
        newdivisor = 0.0
        summand = 0.0
        for A in xrange (self.S, self.J):
            lsummand = self.factor + log(x[0]) * self.S - \
                       lpochham (x[1], int (self.J)) + K [int (A)] + \
                       A * log (x[1]) - lpochham (x[0], int (A)) - divisor
            ...

    def get_kda (self):

        unq_abd = sorted (list (set (self.abund))) # g
        len_unq = len (unq_abd)                    # NDA = i 
        frq_unq = table (self.abund)               # f
        max_a   = max(self.abund)                  # MaxA
        phi = [0] * (max_a +1)
        for s in xrange (self.S):
            phi [self.abund[s]] += 1
        T = [[] for _ in xrange (len_unq)]
        T[0] = [0] * (unq_abd[0] + 1)
        T[0][0] = 0
        T[0][1] = 1
        if unq_abd [0] != 1:
            ls2 = [0] * (unq_abd[0] + 1)
            ls2[1] = 1
            for n in xrange (2, unq_abd [0]+1):
                ls1 = [0] * (n + 1)
                for im in xrange (n):
                    ls1 [im] = ls2 [im]
                ls1 [n] = 0
                for im in xrange (2, n+1):
                    ls2[im] = ls1[im] + float (ls1[im-1] * (im - 1)) / (n - 1)
            for im in xrange (2, unq_abd[0]+1):
                T[0][im] = ls2[im]

        for _in in xrange (1, len_unq):
            T[_in] = [0] * (unq_abd [_in] + 1)
            T[_in] [1] = 1
            ls2 = [0] * (unq_abd[_in] + 1)
            for im in xrange (unq_abd [_in-1]+1):
                ls2[im] = T[_in-1][im]
            for n in xrange (unq_abd[_in-1]+1, unq_abd[_in]+1):
                ls1 = [0] * (n + 1)
                for im in xrange (n):
                    ls1 [im] = ls2 [im]
                ls1 [n] = 0
                for im in xrange (2, n + 1):
                    ls2 [im] = ls1 [im] + float (ls1 [im - 1] * (im - 1)) / (n - 1)
            for im in xrange (2, unq_abd [_in] + 1):
                T [_in][im] = ls2 [im]

        K = [0] * self.J
        poly2 = K[:]
        K [0] = 1
        degree  = 0
        for i in  xrange (len_unq):
            for j in xrange (frq_unq[i]):
                for nn in xrange (degree+1):
                    for mm in xrange (1, unq_abd [i]+1):
                        if K[nn] > 0:
                            poly2 [nn+mm] += T[i][mm]*K[nn]
                degree += unq_abd[i]
                for nn in xrange (degree+1):
                    K [nn] = float (poly2[nn])/10**(4500.0/self.S)
                    poly2[nn] = 0.0

from subprocess import Popen
from os import listdir
from cPickle import load


def parthy_cmd (parthy_path, abd):
    '''
    run parthy just has function name!
    '''
    return '%s -i %s -s 1 -n 0 -N 100000 -J %i -t1 %f -I1 %f -L1 %1.2f\n' % \
           (parthy_path, abd.data_path[:-4],
            abd.J_tot, abd.theta, float(abd.J)/10, abd.lnl)

def main():
    """
    main function
    infile = '/home/francisco/tools/parthy/bci.txt'
    """
    dataset = 'my_mask'
    PATH = '/home/francisco/project/repetitiones/dataset/pickle_%s/' % (dataset)
    parthy_path = '/home/francisco/bin/parthy'
    parthy_dir = '/home/francisco/project/repetitiones/dataset/parthy_ewens/'
    
    # repeat dict
    repeats    = {}
    # dict randomization
    chr_len = {}
    for spe in listdir (PATH):
        if spe.startswith (('Ornithorhynchus', 'Phaeodactylum')):
            continue
        if not spe.endswith ('_repeats.pik'):
            continue
        spe = '_'.join(spe.split('_')[:-1])
        try:
            print spe
            repeats.update (load (open (PATH + spe + '_repeats.pik')))
            chr_len.update (load (open (PATH + spe + '_length.pik' )))
        except IOError:
            continue

    print 'Computing abundances for:'
    abundances = []
    Popen (['mkdir', '-p', parthy_dir]).wait()
    for spe in repeats:
        print '    '+spe
        Popen (['mkdir', '-p', parthy_dir + spe.replace (' ','_')]).wait()
        for crm in chr_len [spe]:
            Popen (['mkdir', '-p',
                    parthy_dir + spe.replace (' ','_') + '/' + crm]).wait()
            infile = parthy_dir + spe.replace (' ',
                                               '_') + '/' + crm + '/reps.txt'
            reps = []
            for rep in repeats[spe]:
                if not crm in repeats[spe][rep] or rep == 'total':
                    continue
                reps.append(repeats[spe][rep][crm])
            open (infile, 'w').write ('\n'.join (map (str, reps)))
            abundances.append ( Abundance (infile,
                                           J_tot=sum(repeats [spe]['total'].values())))

    cola = open ('Parthy_ewens.q', 'w')
    for abd in abundances:
        cola.write (parthy_cmd (parthy_path, abd))
    cola.close()

    

if __name__ == "__main__":
    exit(main())

