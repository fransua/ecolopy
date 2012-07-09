#!/usr/bin/python
"""
22 Sep 2011


"""

__author__  = "Francois-Jose Serra"
__email__   = "francois@barrabin.org"
__licence__ = "GPLv3"
__version__ = "0.0"


from time import sleep, time


print """
When starting with ECOLOGY package the first thing is to load it:

  from ecolopy_dev import Community

"""

from ecolopy_dev import Community

sleep (1)

print """
Community is a class, and derived objects represent simply a distribution of
species abundance , with associated function in order to calculate descriptive
statistics or to fit it to evolutionary models.

The data needed to create this object consists on a list of values corresponding
to the abundance of each species. We can either give to Community a python list
of values:
  Community ([1, 4, 4, 12, 54, 68, 32, 15])
or the path to a file containing those values:
  Community ('whatever_path/mydata.txt')

mydata.txt would contain the same list of values, one per row:
  1
  4
  4
  12
  54
  68
  32
  15

... type 'enter' to continue
"""

raw_input()

print """
In the next step we are going to create an object 'Community' that will
represent the distribution of abundance of the well known/studied BCI dataset.

We are going to load this object under the name 'com':

>>> com = Community ('bci_full.txt')

"""

com = Community ('../dataset_trial/bci_full.txt')

print """
In order to see quickly how does this abundance looks like, we can use the print
command:

>>> print com
  
"""

print com

print """
... type 'enter' to continue
"""
raw_input()

print """

 * Number of individuals correspond to the total of the given community
 * Number of species should correspond to the number of element in your input
   list, or to the number of line in your input file
 * Shannon entropy is computed according to:
    ____
    \\
     \\   p(Xi) * log (pXi)
     /
    /___
   i=1 -> n

   X being the number of individuals for each species
   n the number of species

 * Metacommunity size: correspond to 3 times the community size if not defined
   by user, we could have write to fix it at a given value instead of default:
      com = Community ('bci_full.txt', j_tot=10000000)
 * Models computed: Community can be associated to an Ecological model, the user
   need first to compute them.
 * Current Model: once computed, we can associate our abundance to a given model
 * theta: given by the model
 * I: given by the model
 * m: given by the model
    
"""


print """
... type 'enter' to continue
"""
raw_input()

print """
Once our distribution of abundances loaded into an Community object, we can try
to fit it into an ecological model like Ewens model that assumes that:
  I  = m / (1 - m) * (J - 1)

we just have to type:

>>> com.fit_model(name='ewens')

this step is usually very fast.
  
"""

com.fit_model(name='ewens')


print """
... type 'enter' to continue
"""
raw_input()

print """
We should now see this model into the "Models computed" when doing:
>>> print com
  
"""

print com


print """
... type 'enter' to continue
"""
raw_input()

print """
to load this model as our current model, just type:
>>> com.set_current_model('ewens')
  
"""

com.set_current_model('ewens')

print """
and print again:
"""

print com


print """
We can see here the estimations of theta m, and I according to Ewens model
"""

print """
... type 'enter' to continue
"""
raw_input()


print """
Once loaded some functions are available, for example we can generate a random
neutral distribution of abundance according to this model:
  com.generate_random_neutral_distribution()
"""

print com.generate_random_neutral_distribution()

print """
EcoloPy package use GMP library in order to deal with huge number, usually we
want to get 'normal' numbers in order to compute mean standard deviation...
using common python packages.
Those numbers are quite ugly but easy convert those into standard integers or
floats:

>>> [int (i)for i in com.generate_random_neutral_distribution()]

"""

print [int (i)for i in com.generate_random_neutral_distribution()]


print """
... type 'enter' to continue
"""
raw_input()
model = com.get_model('ewens')
print """
Now we have fit our abundance to Ewens model, a summary of corresponding
parameters are available through the print function, but each of them can also
be reach like this:
>>> com.theta
%s
>>> com.m
%s
>>> com.I
%s

the likelihood of a model can be retrieved like this:
>>> model = com.get_model('ewens')
>>> model.likelihood(com,theta)
%s

or

using the model:

>>> model = com.get_model('ewens')
>>> model.lnL
%s

""" % (com.theta, com.m, com.I, model.likelihood(com.theta),
       com.get_model('ewens').lnL)

t0 = time()
print """
Now we can run an other model like the one proposed by Etienne (2005), just type
>>> com.fit_model(name='etienne')

... it will take around 3-5 minutes (verbosity can be avoided like this:
                                     com.fit_model(name='etienne'))
"""

com.fit_model(name='etienne')
took = time()-t0
print '   ...it took:', int(took/60), 'min', int (took%60), 'sec'


print """
... type 'enter' to continue
"""
raw_input()

com.set_current_model ('etienne')

print """
we have run now one new model, Community object are able to decide though LRT
which has a better likelihood (this is done with chi square test with one degree
of freedom, corresponding to the estimation of parameter m):

>>> com.lrt ('ewens', 'etienne')
%s

this means that etienne model has a significantly better fit than ewens.
(its likelihood is: com.get_model('etienne').lnL -> %s )

Now we know it is the best model we can use it for our dataset:
>>> com.set_current_model ('etienne')
and see the estimations of its parameters:
>>> print com
%s

""" % (com.lrt ('ewens', 'etienne'), com.get_model('etienne').lnL,
       com.__str__())



print """
... type 'enter' to continue
"""
raw_input()


print """
Once chosen and set as current model our best-fit model, we can compute
random distributions of abundance in order to compare them with our real data.
Two examples of random distributions:
>>> com.generate_random_neutral_distribution()
%s
>>> com.generate_random_neutral_distribution()
%s
""" % (com.generate_random_neutral_distribution(),
       com.generate_random_neutral_distribution())

t0 = time()
print """
This is exactly the starting point in order to compute a neutrality test.
First we should generate many random distributions, and compute their shannon
entropy. To finally compare the distribution of shannon entropy to our real
data.

All this is done with the test_neutrality function,  by default only
100 generations will be computed, so it is a good thing to change default
options, this function will return the p-value of rejection of neutrality:
>>> com.test_neutrality(gens=1000)
%s

(still this is a quick example, number of generation around 100000
should be used)
""" % (com.test_neutrality(gens=1000))
took = time()-t0
print '   ...it took:', int(took/60), 'min', int (took%60), 'sec'

print """
... type 'enter' to continue
"""
raw_input()


print """
This is it for the main function of the package, however user will find nice.
ability to store an abundance object:
>>> com.dump_abundance('stored_bci.pik')
to load it afterward, (it will avoid from running again computation of K(D,A))
>>> com = Community ('stored_bci.pik')


....

The End.
"""


