#!/usr/bin/python
"""
04 Jul 2011


"""

__author__  = "Francois-Jose Serra"
__email__   = "francois@barrabin.org"
__licence__ = "GPLv3"
__version__ = "0.13"

from scipy.stats    import chisqprob#, lognorm
from math           import log, exp
from cPickle        import dump, load
from os.path        import isfile, exists
from sys            import stdout, stderr
from numpy          import arange, mean, std
from warnings       import warn

from ecolopy_dev.utils  import shannon_entropy
from ecolopy_dev.models import *

class Community (object):
    '''
    Main community class.

    :argument data: containing species count, can be given in one of those format:
    
       * python list, each element being a species count
       * text file containing one species count per line
       * pickle file containing community object

    :argument None j_tot: is the size of the metacommunity, if not defined j_tot = J * 3

    :returns: an Community object

    **Examples:**
    ::

      # python list
      abd = Community([1, 1, 2, 3, 4, 8, 12])
      # text file
      abd = Community("bci.txt", j_tot=256987)
      # finally a saved Community object
      abd = Community("abundance_bci_first_try.pik")
    
    '''

    def __init__ (self, data, j_tot=None):
        """
        creation of Community object
        """
        self.__models  = {}
        if type (data) != list:
            self.data_path = data
            data = self._parse_infile ()
        if not data is None:
            self.abund     = [x for x in sorted (data [:])]
        self.J         = sum (self.abund)
        self.S         = len (self.abund)
        self.j_tot     = j_tot if j_tot else self.J * 3
        self.shannon   = shannon_entropy (self.abund, self.J)
        self.__current_model = None
        
    def __str__(self):
        """
        for print
        """
        return '''Community (object)
    Number of individuals (J) : %d
    Number of species (S)     : %d
    Shannon's index (shannon) : %.4f
    Metacommunity size (j_tot): %d
    Models computed           : %s
    Model loaded              : %s
    ''' % (self.J, self.S, self.shannon, self.j_tot,
           ', '.join (self.__models.keys()),
           (self.__current_model.__class__.__name__ + \
            '\n' + ' ' * 8 + ('\n' + ' ' * 8).join(
                [l for l in str(self.__current_model).split('\n')[2:]]
            )) if self.__current_model else "None")


    def rsa_ascii (self, width=100, height=50, pch='o'):
        """
        Draw Relative Species Abundances curve (ASCII format).

        :argument 100 width: width in term of characters
        :argument 100 height: height in term of characters
        :argument o pch: dot character

        :returns: string corresponding to plot
        
        """
        dots = sorted([float(x) for x in self.abund],reverse=True)
        S = float(self.S)
        J = float(self.J)
        rabd = []
        for d in dots:
            rabd.append (log (100*d/J))
        y_arange = sorted(arange(min(rabd), max(rabd), abs (min(rabd)- max(rabd))/height), reverse=True)
        x_arange = sorted(arange(0, S, S/width))
        y = 0
        x = 0
        spaces = ''
        graph = '\n(%) Relative\nAbundances\n\n'
        graph += '{0:<7.4f}+'.format (exp(max(rabd)))
        for i, d in enumerate (rabd):
            if d < y_arange[y]:
                graph += '\n'
                if not (y)%5 and y != 0:
                    graph += '{0:<7.4f}+'.format (exp(y_arange[y]))
                else:
                    graph += ' '*7 + '|'
                graph += spaces
                while d < y_arange[y]:
                    y+=1
            if i > x_arange[x]:
                graph += pch
                spaces += ' '
                if  x+1 >= width:
                    break
                x+= 1
        graph += '\n'
        graph += ' 1/inf ' + ''.join(
            ['+' if not x%5 else '-' for x in xrange(width+1)]) + '\n'
        graph += ' '*7 + ''.join(
            ['{0:<5}'.format(int(x_arange[x])) for x in xrange(0,width,5)]
        ) + str( int(x_arange[-1])) + '\n\n'
        graph += ' '*7 + '{0:^{1}}'.format('Species rank', width)
        return graph

    def draw_rsa(self, outfile=None, filetype=None, size=(15,15)):
        """
        Draw Relative Species Abundances curve.

        :argument None outfile: path were image will be saved, if none, plot will be shown using matplotlib GUI
        :argument None filetype: pdf or png
        :argument (15,15) size: size in inches of the drawing

        """
        import pylab
        try:
            import matplotlib.ticker
        except ImportError:
            warn("WARNING: matplotlib not found, try rsa_ascii function instead")
            
        pylab.grid(alpha=0.4)
        y  = []
        x  = []
        for rank,  abd in enumerate(sorted(self.abund, reverse=True,
                                           key=float)):
            x.append (rank+1)
            y.append (100*float (abd) / float(self.J))
        if len (y) == 1:
            raise Exception ('list of abundances is too short man')
        mark = 'o'
        lcol = 'grey'
        mcol = 'black'
        lsty = '-'
        pylab.plot(x, y, marker=mark, ms=5, color=lcol, lw=3, mfc=mcol, ls=lsty)
        # title legend...
        maxX = len (x)
        pylab.yscale('log')
        pylab.xticks(range (0,maxX+1,5), range (0,maxX+1,5), rotation=90)
        pylab.title('Ranked Species Abundance')
        pylab.xlabel('Species rank number (rank according to abundance)')
        pylab.ylabel('Relative abundance percentage of each species')
        pylab.xlim(1, len(x)+1)
        pylab.ylim(log(1.0/float(self.J)), max(y)*1.5)
        ax = pylab.gca()
        ax.yaxis.set_major_formatter(matplotlib.ticker.FormatStrFormatter('%.4f'))
        F = pylab.gcf()
        DPI = F.get_dpi()
        F.set_size_inches (size)
        if outfile:
            if filetype == 'pdf':
                F.savefig (outfile, dpi=DPI+30, filetype='pdf')
            elif filetype == 'png':
                F.savefig (outfile, dpi=DPI+30, filetype='png')
            pylab.close()
        else:
            pylab.show()
        

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

    def fit_model (self, name="ewens", **kwargs):
        """
        Fit Community to model.
        Extra arguments can be pssed depending on the model to fit. See doc of its corresponding optimize function.

        :argument Ewens name: name of the model, between Etienne, Ewens or Log-normal

        :return: model
        """
        if name == 'etienne':
            self.__models['etienne']   = EtienneModel(self, **kwargs)
        elif name == 'ewens':
            self.__models['ewens']     = EwensModel(self, **kwargs)
        elif name == 'lognormal':
            self.__models['lognormal'] = LognormalModel(self, **kwargs)

    def set_model (self, model):
        """
        add one model computed externally to the computed models of current Community object

        :argument model: model object
        
        """
        name  = model.__class__.__name__.replace('Model', '').lower()
        self.__models[name] = model


    def get_model (self, name):
        """
        :argument name: name of a computed model
        :returns: a EcologicalModel object corresponding to one of the computed models
        """
        if name in self.__models:
            return self.__models[name]
        else:
            return None


    def set_current_model (self, name):
        """
        set on model as default/current model.

        :argument name: model name of precomputed model
        
        """
        self.__current_model = self.get_model(name)
        for key in ['theta', 'I', 'm']:
            try:
                self.__dict__[key] = self.__current_model._parameters[key]
            except KeyError:
                self.__dict__[key] = None
        self.__dict__['lnL'] = self.__current_model.lnL
                

    def test_neutrality (self, model='ewens', gens=100, full=False, fix_s=False,
                         tries=1000, method='shannon', verbose=False):
        '''
        test for neutrality comparing Shannon entropy
        if (Hobs > Hrand_neut) then evenness of observed data is
        higher then neutral

        :argument ewens model: model name otherwise, current model is used
        :argument 100 gens: number of random neutral distributions to generate
        :argument False full: also return list of compared values (H or lnL) for simulated abundances
        :argument False fix_s: decide whether to fix or not the number of species for the generation of random neutral communities
        :argument False tries: in case S is fixed, determines the number of tries in order to obtain the exact same number of species as original community. In case The number of tries is exceeded, an ERROR message is displayed, and p-value returned is 1.0.
        :argument shannon method: can be either "Shannon" for comparing Shannon's entropies or "loglike" for comparing log-likelihood values (Etienne 2007). Last method is much more computationally expensive, as likelihood of neutral distribution must be calculated.
        :returns: p_value if full=True also returns Shannon entropy (or likelihoods if method="loglike") of all random neutral abundances generated
        '''
        fast_shannon = lambda abund: sum ([-spe * log(spe) for spe in abund])
        pval = 0
        inds = self.S if 'LognormalModel' in repr(model) else self.J
        model = self.get_model (model)
        if not model:
            warn("WARNING: Model '%s' not found.\n" % model)
            return None
        neut_h = []
        for _ in xrange (gens):
            if verbose:
                stdout.write ("\r  Generating random neutral abundances %s out of %s" \
                              % (_+1, gens))
                stdout.flush ()
            if fix_s:
                for _ in xrange (tries):
                    tmp = model.random_community(inds)
                    if len(tmp) == self.S:
                        break
                else:
                    stderr.write('ERROR: Unable to obtain S by simulation')
                    if full:
                        return float('nan'), []
                    return float('nan')
            else:
                tmp = model.random_community(inds)
            l_tmp = sum (tmp)
            if method == 'shannon':
                neut_h.append ((fast_shannon (tmp) + l_tmp*log(l_tmp))/l_tmp)
                pval += neut_h[-1] < self.shannon
            elif method == 'loglike':
                tmp = Community(tmp)
                tmp._get_kda(verbose=False)
                neut_h.append (tmp.etienne_likelihood([model.theta,model.m]))
                pval += neut_h[-1] < model.lnL
        if verbose:
            stdout.write ('\n')
        if full:
            return float (pval)/gens, neut_h
        return float (pval)/gens
    

    def _parse_infile (self):
        '''
        TODO: use other columns
        parse infile and return list of abundances
        infile can contain several columns
        '''
        abundances = []
        if exists(self.data_path):
            lines = open (self.data_path).readlines()
        else:
            lines = self.data_path.split('\n')
        if lines[0].strip() == '(dp1':
            self.load_community (self.data_path)
            return None
        for line in lines:
            if line:
                abundances.append (int (line.strip().split('\t')[0]))
        return abundances


    def dump_community (self, outfile, force=False):
        '''
        save params and kda with pickle
        force option is for writing pickle with no consideration
        if existing

        :argument outfile: path of the outfile
        :argument False force: overwrite existing outfile
        
        '''
        if isfile (outfile) and not force:
            self.load_community (outfile)
        # changed, kda no longer here
        #try:
        #    self.__models['KDA']   = self._kda[:]
        #except TypeError:
        #    self.__models['KDA'] = None
        self.__models['ABUND'] = self.abund[:]
        dump (self.__models, open (outfile, 'w'))
        del (self.__models['ABUND'])


    def load_community (self, infile):
        '''
        load params and kda with pickle from infile
        WARNING: do not overwrite values of params.

        :argument infile: path of the outfile

        '''
        old_params = load (open (infile))
        old_params.update(self.__models)
        self.__models = old_params
        self.abund = self.__models['ABUND']
        del (self.__models['ABUND'])
        #if 'KDA' in self.__models:
        #    self._kda  = self.__models['KDA']
        #    del (self.__models['KDA'])
        #else:
        #    self._kda = None


    def generate_random_neutral_distribution (self, model=None, J=None):
        """
        return distribution of abundance, according to a given model
        and a given community size.
        if none of them are given, values of current Community are used.

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
            J = self.S if 'LognormalModel' in repr(model) else self.J
        return model.random_community (J)

