
.. moduleauthor:: Francois Serra
.. currentmodule:: ecolopy

Just after counting species abundances in an ecosystem
******************************************************

.. contents::

Load Abundance
==========

Abundance is a class, and derived objects represent simply a distribution of
species abundance , with associated function in order to calculate descriptive
statistics or to fit it to evolutionary models.

The data needed to create this object consists on a list of values corresponding
to the abundance of each species. We can either give to Abundance a python list
of values:
::

  Abundance ([1, 4, 4, 12, 54, 68, 32, 15])

or the path to a file containing those values:
::

  Abundance ('whatever_path/mydata.txt')

mydata.txt would contain the same list of values, one per row:

::

  1
  4
  4
  12
  54
  68
  32
  15

In the next step we are going to create an object 'Abundance' that will
represent the distribution of abundance of the well known/studied BCI dataset.

We are going to load this object under the name 'abd':

::

  from ecolopy import Abundance
  abd = Abundance ('bci_full.txt')

In order to see quickly how does this abundance looks like, we can use the print
command:

::

  print abd

  #  Abundance (object)
  #      Number of individuals (J) : 21457
  #      Number of species (S)     : 225
  #      Shannon entropy (shannon) : 4.2704
  #      Metacommunity size (j_tot): 64371
  #      Models computed           : 
  #      Current model             : None
  #      Theta                     : None
  #      I                         : None
  #      m                         : None

With:

* Number of individuals correspond to the total of the given community
* Number of species should correspond to the number of element in your input list, or to the number of line in your input file
* Shannon entropy is computed according to:

.. math::
  :nowrap:

  \begin{eqnarray}
    H(X) = \sum_{0\le i\le n} p(x_i) * log (p(x_i))
  \end{eqnarray}
    X being the number of individuals for each species and n the number of species

* Metacommunity size: correspond to 3 times the community size if not defined by user, we could have write to fix it at a given value instead of default:

::

  abd = Abundance ('bci_full.txt', j_tot=10000000)

* Models computed: Abundance can be associated to an Ecological model, the user
  need first to compute them.
* Current Model: once computed, we can associate our abundance to a given model
* theta: given by the model
* I: given by the model
* m: given by the model


Fit to ecological model
***********************

Once our distribution of abundances loaded into an Abundance object, EcoloPy proposes a set of of ecological models that we can try to fit to our data.

Ewens model
===========

This model assumes that:

.. math::
  :nowrap:

  \begin{eqnarray}
    I  = \frac{m}{(1 - m) * (J - 1)}
  \end{eqnarray}

we just have to type:

::

  abd.ewens_optimal_params()

this step is usually very fast.

to load this model as our current model, just type:

::

  abd.set_current_model('ewens')

  print abd

  #  Abundance (object)
  #      Number of individuals (J) : 21457
  #      Number of species (S)     : 225
  #      Shannon entropy (shannon) : 4.2704
  #      Metacommunity size (j_tot): 64371
  #      Models computed           : ewens
  #      Current model             : ewens
  #      Theta                     : 34.962254203932339
  #      I                         : 17.494565308269266
  #      m                         : 0.00081470508933989701


Etienne model
=============

Now we can run an other model like the one proposed by Etienne (2005), just type

::

  abd.etienne_optimal_params()
  abd.set_current_model ('etienne')
  # Getting K(D,A) according to Etienne 2005 formula:
  #   Getting some stirling numbers...
  #     1000 of 1717.0, size: 3145976
  #   Computing K(D,A) at species 108 out of 108

  print abd
  # Abundance (object)
  #     Number of individuals (J) : 21457
  #     Number of species (S)     : 225
  #     Shannon entropy (shannon) : 4.2704
  #     Metacommunity size (j_tot): 64371
  #     Models computed           : ewens, etienne
  #     Current model             : etienne
  #     Theta                     : 47.6743190606
  #     I                         : 2211.0866293821641
  #     m                         : 0.0934245377983


Lognormal model
===============

In ecolopy is also implemented log normal model:

::

  abd.lognorm_optimal_params()
  abd.set_current_model ('lognorm')

  print abd
  # Abundance (object)
  #     Number of individuals (J) : 21457
  #     Number of species (S)     : 225
  #     Shannon entropy (shannon) : 4.2704
  #     Metacommunity size (j_tot): 64371
  #     Models computed           : lognorm, ewens, etienne
  #     Current model             : lognorm
  #     Theta                     : 3.14269458985
  #     I                         : 1.78719175872
  #     m                         : None


Comparing Models
****************

Browsing parameters
===================

Now we have fit our abundance to some models, a summary of corresponding
parameters are available through the print function, but each of them can also
be reach like this:

::

  abd.theta
  # 34.962254203932339
  abd.m
  # 0.00081470508933989701
  abd.I
  # 17.494565308269266

To get the likelihood, we can use using the model Object linked to our abundance:
::

  model = abd.get_model('ewens')
  model.lnL
  # 318.84864864917472

All those values coorespond to the current model (in this case ewens).


Searching for best model
========================

We have run now several models, within which the nested models ewens and etienne, for those we can run a likelihood ratio test using the lrt function.
This function will compute a chi square test for 2 times the difference in likelihoods, with one degree of freedom (corresponding to the estimation of parameter m):

::

  abd.lrt ('ewens', 'etienne')
  # 6.80784682569e-06


Generate random distribution
****************************

In order to compare our distribution of abundance to the expected one according to a specific model, we can generate random neutral distributions.
By default EcoloPy will use the parameters of the current model but this can be change passing to the function the name of the wanted model:

::

  abd.generate_random_neutral_distribution(model='etienne')
  # [mpfr('17.0'), mpfr('867.0'), mpfr('397.0'), mpfr('184.0'), mpfr('71.0'), 
  # ...
  # ...
  # mpfr('2.0'), mpfr('2.0'), mpfr('1.0'), mpfr('1.0'), mpfr('1.0'), mpfr('1.0')]


Note: EcoloPy package use GMP library in order to deal with huge number, usually we
want to get 'normal' numbers in order to compute mean standard deviation...
using common python packages.
Those numbers are quite ugly but easy convert those into standard integers or
floats:

::

  abd.generate_random_neutral_distribution()
  # [mpfr('17.0'), mpfr('867.0'), mpfr('397.0'), mpfr('184.0'), mpfr('71.0'), 
  # ...
  # ...
  # mpfr('2.0'), mpfr('2.0'), mpfr('1.0'), mpfr('1.0'), mpfr('1.0'), mpfr('1.0')]

  # or in order to get floats:
  [int (i)for i in abd.generate_random_neutral_distribution()]
  # [273, 263, 461, 754, 1140, 163, 67, 113, 1014, 407, 1496, 1395, 405, 534, 1435, 260, 
  # ...
  # ...
  # 3, 2, 1, 1, 3, 1, 1, 1, 1, 1, 1]



Saving/Loading Abundance object
*******************************

Once done EcoloPy allow user to save Abundance object and EcologicalModels object into cPikle with dum_abundance and load_abundance functions.

::

  # save it
  abd.dump_abundance('stored_bci.pik')
  # (re)load it
  abd = Abundance ('stored_bci.pik')
