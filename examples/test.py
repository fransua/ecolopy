

from ecolopy import Abundance
from ecolopy.utils import draw_shannon_distrib

abd = Abundance('test.txt')
abd.lognorm_optimal_params()
abd.set_current_model('lognorm')
pval, neut_h = abd.test_neutrality (model='lognorm',
                                    gens=1000, give_h=True)
draw_shannon_distrib(neut_h,abd.shannon)

abd2 = Abundance('test2.txt')
abd2.lognorm_optimal_params()
abd2.set_current_model('lognorm')
pval, neut_h = abd2.test_neutrality (model='lognorm',
                                    gens=1000, give_h=True)
draw_shannon_distrib(neut_h,abd2.shannon)


bci = Abundance('../dataset_trial/bci_full.txt')
bci.lognorm_optimal_params()
bci.set_current_model('lognorm')

pval, bci_neut_h = bci.test_neutrality (model='lognorm',
                                    gens=1000, give_h=True)


print pval
draw_shannon_distrib(bci_neut_h,bci.shannon)

from scipy.stats import kstest, lognorm
import numpy as np
import pylab as plb

a = abd.abund
b = bci.abund
c = abd2.abund

a = np.array ([float(i) for i in a])
b = np.array ([float(i) for i in b])
c = np.array ([float(i) for i in c])

log_a = np.array ([np.log (i) for i in a])
log_b = np.array ([np.log (i) for i in b])
log_c = np.array ([np.log (i) for i in c])

log_a = (log_a-np.mean (log_a))/np.std (log_a)
log_b = (log_b-np.mean (log_b))/np.std (log_b)
log_c = (log_c-np.mean (log_c))/np.std (log_c)

print kstest (log_a, 'norm')
print kstest (log_b, 'norm')
print kstest (log_c, 'norm')

plb.hist (b)
plb.hist (log_b, bins=20)
plb.hist (a, bins=100)
plb.hist (log_a, bins=10)

shape, loc, scale = lognorm.fit(a)
rnd_a = lognorm.rvs(shape, scale=scale, loc=loc, size=len(a))
plb.hist(rnd_a, bins=20, alpha=0.5)
plb.hist(a, bins=20, color='r', alpha=0.5)

shape, loc, scale = lognorm.fit(c)
rnd_c = lognorm.rvs(shape, scale=scale, loc=loc, size=len(c))
plb.hist(rnd_c, bins=30, alpha=0.5)
plb.hist(c, bins=30, color='r', alpha=0.5)

shape, loc, scale = lognorm.fit(b)
rnd_b = lognorm.rvs(shape, scale=scale, loc=loc, size=len(b))
plb.hist(rnd_b, bins=20, alpha=0.5)
plb.hist(b, bins=20, color='r', alpha=0.5)

mean (b)
shape = std (b)

plb.show()
