

Introduction
************

.. contents::


What is EcoloPy?
================

EcoloPy directly draws inspiration from the UNTB package in R [Hankin2007]_, to infer neutrality of an ecosystem distribution of species. It is able to fit sampling data of a community into a neutral model and implements functions to interact and represent the original data and the results of statistical tests.
Ewens and Etienne neutral models are implemented.

Why an other package?
=====================

Several packages or programs were already developed in order to deal with species abundances data, implementing statistical functions in order to fit data in ecological models and even able to test for neutrality [JabotChave2011]_ [Etienne2007]_ [Hankin2007]_. However none of those programs were able to deal with genomic data, with abundances reaching the million of individuals. The main computational bottleneck, or in this case cul-de-sac, being the calculation of Etienne’s sampling formula [Etienne2005]_ where computation of K(D, A) uses stirling numbers.

In order to adapt the algorithm already implemented, to genomic dataset we developed the Ecolopy package, that, as a main point, uses the GMP [Granlund2000]_ and MPFR [Fousse2007]_ libraries through GMPY biding [Martelli2007]_. Other improvements specific to Ecolopy, and needed for dealing with genomic dataset where implemented, mainly:
  * usage of recurrence function of stirling numbers, and building of table of stirling numbers (already used by [JabotChave2011]_)
  * table of stirling numbers is reduced on the fly keeping only numbers needed, in order to save memory.
  * model optimization using different optimization strategies from Scipy [Jones2001]_.

In top of those necessary technical improvements, EcoloPy presents one main advantage as it is entirely written in Python [vanRossumdeBoer1991]_, a programming language that offers a strong support for integration with other languages and tools, and whose popularity is raising among the bioinformatics community [Bassi2007]_. Ecolopy is still a fully ripened package as the number of functions and ecological models is still low. But the package was designed to provide a scalable program architecture.

The ecological models currently implemented are:
  * Hubbell's implementation of Ewens sampling formula [Hubbell2001]_.
  * Etienne's sampling formula [Etienne2005]_.
  * Log-normal.



References
==========

.. [Bassi2007] S Bassi, A primer on python for life science researchers. PLoS computational biology 3(11) (2007), p. e199.

.. [Etienne2005] Rampal S Etienne, A new sampling formula for neutral biodiversity. Ecology Letters 8(3) (2005), 253–260.

.. [Etienne2007] Rampal S Etienne, A neutral sampling formula for multiple samples and an 'exact' test of neutrality. Ecology letters 10(7) (2007), 608–18.

.. [Fousse2007] Laurent Fousse, Guillaume Hanrot, Vincent Lefevre, Patrick Pelissier, and Paul Zimmermann, MPFR: A multiple-precision binary floating-point library with correct rounding. ACM Transactions on Mathematical Software (TOMS) 33(2) (2007), p. 13.

.. [Granlund2000] Torbjorn Granlund, GMP: The GNU Multiple Precision Arithmetic Library, 2000.

.. [Hankin2007] R.K.S. Hankin, Introducing UNTB, an R package for simulating ecological drift under the unified neutral theory of biodiversity. Journal of Statistical Software 22(12) (2007), 1–15.

.. [Hubbell2001] Stephen P Hubbell, The unified Neutral Theory of Biodiversity and Biogeography. Princeton University Press, 2001.

.. [JabotChave2011] Franck Jabot and Jerome Chave, Analyzing Tropical Forest Tree Species Abundance Distributions Using a Nonneutral Model and through Approximate Bayesian Inference. The American naturalist 178(2) (2011), E37–47.

.. [Jones2001] Eric Jones, Travis Oliphant, Pearu Peterson, and Others, Scipy: Open source scientific tools for Python, 2001.

.. [Martelli2007] Alex Martelli, GMPY Multiprecision arithmetic for Python, 2007.

.. [vanRossumdeBoer1991] G van Rossum and J de Boer, Interactively testing remote servers using the python programming language. CWI Quarterly 4(4) (1991), 283–303.
