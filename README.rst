=======
EcoloPy
=======

**a Python package to test for Unified Neutral Theory of Biodiversity and biogeography in Large dataset.**
EcoloPy deal with species abundance and diversity data. It allows to fit this data into ecological models and implements some statistical tests.
It is able to deal with large dataset like genetic element content of mammalian genomes.

The package was implemented in the context of a study of the distribution of genetic elements in Eukaryotic genomes `[Serra2013]`_.


Documentation
=============

http://fransua.github.io/ecolopy/


.. _[Serra2013]:

Citing EcoloPy
--------------

Please refer to the following article if you use this package to analyze or process results that are part of a published work:

Serra, F., Becher, V., & Dopazo, H. (2013). 
**Neutral Theory Predicts the Relative Abundance and Diversity of Genetic Elements in a Broad Array of Eukaryotic Genomes.** 
PloS One, 8(6), e63915. doi:`10.1371/journal.pone.0063915 <http://dx.plos.org/10.1371/journal.pone.0063915>`_


PREREQUISITE
============

PYTHON
------

This is to install python packages easily.

::

  sudo apt-get install python-setuptools


MPFR 
-----

http://www.mpfr.org/

::

  cd /tmp/
  wget http://www.mpfr.org/mpfr-current/mpfr-3.1.2.tar.gz
  tar xzvf mpfr-3.1.2.tar.gz
  cd mpfr-3.1.2/
  ./configure
  make
  make check
  sudo make install

GMP
---

http://gmplib.org/

::

  cd /tmp/
  wget ftp://ftp.gmplib.org/pub/gmp/gmp-5.1.2.tar.bz2
  tar xjvf gmp-5.1.2.tar.bz2
  cd gmp-5.1.2/
  ./configure
  make
  make check
  sudo make install

MPC
---

http://www.multiprecision.org/

::

  cd /tmp/
  wget http://www.multiprecision.org/mpc/download/mpc-1.0.1.tar.gz
  tar xzvf mpc-1.0.1.tar.gz
  cd mpc-1.0.1/
  ./configure
  make
  make check
  sudo make install


GMPY2 
-----

for huge numbers.
::

  sudo apt-get install python2.7-dev libmpfr-dev libgmp3-dev libmpc-dev
  sudo easy_install gmpy2


SCIPY
-----

for optimization of log likelihood

::

  sudo easy_install scipy


NUMPY
-----

::

  sudo easy_install numpy


MATPLOTLIB
----------

Not compulsory, only to draw graphs

::

  sudo easy_install matplotlib


INSTALL:
========

::

  sudo python setup.py install

test usage:

::

  python test_all.py full


TODO:
-----
* do something with metacommunity


Travis build
------------

*Broken* -> problem with scipy install at travis.

.. |build-status|
   image:: https://secure.travis-ci.org/fransua/ecolopy.png
           ?branch=master
   :target: http://travis-ci.org/fransua/ecolopy
   :alt: Build Status
