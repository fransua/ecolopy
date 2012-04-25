#!/usr/bin/python
"""
27 Dec 2011


"""

from optparse import OptionParser
from ecolopy_dev import Abundance
import pylab
import numpy as np
from scipy.interpolate import spline
from ecolopy_dev.utils import fast_etienne_likelihood

__title__ = 'Ecolopy-UNTBGen'

def draw_rsa(abund, out, smooth=True, labels=False, pdf=False):
    pylab.grid(alpha=0.4)
    maxX = 0
    y  = []
    x  = []
    total = sum (abund.values())
    for spe in abund:
        x.append ([abund[spe], spe])
        y.append (float (x[-1][0]) / total)
    y = np.array (sorted (y, reverse=True))
    #x = [x[1] for x in sorted (x, key = lambda x: x[0], reverse=True)]
    if len (y) == 1:
        raise Exception ('list of abundances is too short man')
    m = 'o'
    c = 'red'
    l = '-'
    if smooth:
        xnew = np.linspace(0,len(x),300)
        y_smooth = spline(range(len (x)),[np.log (i) for i in y],
                          xnew, order=3, kind='smoothest')
        prev_i = y_smooth[0]
        for j, i in enumerate(y_smooth[:]):
            if prev_i < i:
                y_smooth[j] = prev_i
            else:
                prev_i = i
    else:
        xnew = range (len (x))
        y_smooth = [np.log (i) for i in y]
    pylab.plot (xnew,  y_smooth,#[log (i) for i in y],
                marker=m, ms=2, lw=2,ls= l, color=c)
    if len (x)> maxX:
        maxX = len (x)
    # title legend...
    if labels:
        pylab.xticks (range(1,len(x)+1), sorted(abund.keys(),key=lambda x: abund[x], reverse=True), rotation=90)
    else:
        pylab.xticks (range (0,maxX+1,5), range (0,maxX+1,5), rotation=90)
    pylab.yticks ([(np.log (float ('1e'+str(i))) if j==0 else np.log(float ('1e'+str(i))) + np.log(float(j)/10)) for i in range (-6,1) for j in range (0,10)],
                  [((float ('1e'+str(i))*100) if j==0 else '') for i in range (-6,1) for j in range (0,10)])
    pylab.title ('Ranked abundance.')
    pylab.xlabel ('Elements number ranked by size')
    pylab.ylabel ('Log percentage of representation of each species')
    pylab.ylim(ymin=min(np.log(y))-1)
    F = pylab.gcf()
    DPI = F.get_dpi()
    F.set_size_inches ( (15.0,17.0))
    if pdf:
        F.savefig (out, dpi=DPI+30, format='pdf')
    else:
        F.savefig (out, dpi=DPI+30, format='png')
    pylab.close()


def draw_shannon_distrib(neut_h, obs_h, out, pval=1):
    '''
    draws distribution of Shannon values for random neutral

    :argument neut_h: list of Shannon entropies corresponding to simulation under neutral model
    :argument obs_h: Shannon entropy of observed distribution of abundance
    
    '''
    neut_h = np.array ([float (x) for x in neut_h])
    obs_h = float (obs_h)
    pylab.hist(neut_h, 40, color='green', histtype='bar', fill=True)
    pylab.axvline(float(obs_h), 0, color='r', linestyle='dashed')
    pylab.axvspan (float(obs_h) - neut_h.std(), float(obs_h) + neut_h.std(),
                    facecolor='orange', alpha=0.3)
    pylab.xlabel('Shannon entropy (H)')
    pylab.ylabel('Number of observations over %s simulations' % (len (neut_h)))
    pylab.title("Histogram of entropies from %s simulations compared to \nobserved entropy (red), deviation computed from simulation, p-val=%.5f" % (len (neut_h),pval))
    F = pylab.gcf()
    DPI = F.get_dpi()
    F.set_size_inches ( (15.0,17.0))
    F.savefig (out, dpi=DPI+30, format='png')
    pylab.close()


def draw_contour_likelihood (abd, out, model=None, theta_range=None,
                             m_range=None, num_dots=100, local_optima=True,
                             write_lnl=False,verbose=True):
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
            if verbose:
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
            if verbose:
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
            pylab.plot(x[j], y[i], color='w', marker='x')

    # define number of color categories... perhaps not the best way
    perc10 = np.percentile(z, 10)
    perc90 = np.percentile(z, 90)
    levels = list (np.arange (perc90, perc10, -(perc90-perc10)/num_dots))[::-1] + \
             list (np.arange (perc90, z.max(), (z.max()-perc90)/num_dots))[1:] + [z.max()]
    pylab.plot(abd.theta, abd.m, color='w', marker='*')
    pylab.contourf (x, y, z, levels)
    pylab.colorbar (format='%.3f')
    pylab.vlines(abd.theta, 0, abd.m, color='w', linestyle='dashed')
    pylab.hlines(abd.m, 0, abd.theta, color='w', linestyle='dashed')    
    if write_lnl:
        fig = pylab.contour (x, y, z, levels)
        pylab.clabel (fig, inline=1, fontsize=10, colors='black', manual=True)
    pylab.title ("Log likelihood of abundance under Etienne model")
    pylab.xlabel ("theta")
    pylab.ylabel ("m")
    pylab.axis (theta_range + m_range)
    F = pylab.gcf()
    DPI = F.get_dpi()
    F.set_size_inches ( (15.0,17.0))
    F.savefig (out, dpi=DPI+30, format='png')
    pylab.close()



def main():
    """
    main function
    """
    opts = get_options()

    abd = Abundance(opts.abund)

    out = open (opts.out + 'result.txt', 'w')

    out.write('\n* ' + str(abd))
    
    if opts.sad:
        abund = {}
        name = 'Species #{0}'
        for i, a in enumerate(sorted(abd.abund, reverse=True)):
            abund[name.format(i+1)] = float(a)
        draw_rsa(abund, opts.out + 'sad.png', smooth=False, labels=False)
        out.write('\n* SAD drawn at sad.png\n')

    best = None
    if 'ewens' in opts.models:
        abd.ewens_optimal_params()
        best = 'ewens'
        out.write ('\n* Ewens model computed:\n')
        abd.set_current_model('ewens')
        out.write('   theta: {0}\n   lnL: {1}\n'.format(abd.theta, abd.lnL))
    if 'etienne' in opts.models:
        abd.etienne_optimal_params()
        best = 'etienne'
        abd.set_current_model('etienne')
        out.write ('\n* Etienne model computed:\n')
        out.write('   theta: {0}\n   m: {1}\n   lnL: {2}\n'.format(abd.theta, abd.m, abd.lnL))

    if opts.lrt:
        if not 'etienne' in opts.models or not 'ewens' in opts.models:
            raise Exception ('can only compute LRT if both etienne and ewens models are computed')

    if opts.lrt:
        pval = abd.lrt('ewens', 'etienne')
        best = 'etienne' if pval<0.05 else 'ewens'
        out.write ('\n* LRT:\n   Best model: {0}, p-val of LRT: {1}\n'.format(best, pval))
        

    if not best:
        raise Exception ('can not perform neutrality test with no model computed')

    abd.set_current_model(best)
    if opts.test:
        lrt, neut_h = abd.test_neutrality(model=best, gens=int (opts.gens),
                                          full=True, method=opts.meth, fix_s=opts.fix_s)
        draw_shannon_distrib(neut_h, abd.shannon if opts.meth=='shannon' else abd.lnL,
                             opts.out + 'test.png', pval=lrt)
        out.write ('\n* Neutrality test under {0} model\n'.format(best))
        out.write ('  (with S {0} and based on comparison of {1} values)\n'.format('fixed' if opts.fix_s else 'free',
                                                                                    'Shannon' if opts.meth=='shannon' else 'Log-likelihood'))
        out.write ('   p-val: {1}\n   figure: test.png\n'.format(best, lrt, ))

    if opts.contour and not 'etienne' in opts.models:
        raise Exception ('Contour only available with etienne model')

    if opts.contour:
        draw_contour_likelihood(abd, opts.out + 'contour.png', model='etienne')
        out.write ('\n* Contour plot drawn at: contour.png\n\n')

    out.close()


def get_options():
    '''
    parse option from command line call
    '''
    parser = OptionParser(
        version=__title__,
        usage="%prog [options] file [options [file ...]]")
    parser.add_option('-i', dest='abund', metavar="PATH",
                      help='path to abundance file')
    parser.add_option('-o', dest='out', metavar="PATH",
                      help='out directory')
    parser.add_option('-g', dest='gens', metavar="INT",
                      default=1000,
                      help='number of simulations neutrality test')
    parser.add_option('-m', dest='models', metavar='LIST',
                      help='''Ecological models.                            
                      e.g.: -m "ewens,etienne"''')
    parser.add_option('--method', dest='meth', default='shannon',
                      help='''kind of neutrality test, must be shannon or loglike''')
    parser.add_option('--lrt', dest='lrt', action='store_true',
                      help='compute likelihood ratio test between models')
    parser.add_option('--test', dest='test', action='store_true',
                      help='compute neutrality test')
    parser.add_option('--sad', dest='sad', action='store_true',\
                      default=False,\
                      help='draw species abundance/diversity curve')
    parser.add_option('--fixs', dest='fix_s', action='store_true',\
                      default=False,\
                      help='fix number of species during simulations')
    parser.add_option('--contour', dest='contour', action='store_true',\
                      default=False,\
                      help='draw contour plot for etienne optimization')
   
    opts = parser.parse_args()[0]
    if not opts.abund or not opts.out:
        parser.print_help()
        exit()
    if not opts.out.endswith('/'):
        opts.out += '/'
    opts.models = opts.models.split(',')
    return opts
        

if __name__ == "__main__":
    exit(main())
