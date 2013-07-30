=======
EcoloPy
=======

PREREQUISITE
============

PYTHON
------

This is to install python packages easily.

::

  sudo apt-get install python-setuptools


MPFR 
----

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
* revise all doc.
* do something with metacommunity

Links
=====

`Travis CI <https://travis-ci.org/#!/tkf/emacs-jedi>`_ |build-status|

.. |build-status|
   image:: https://secure.travis-ci.org/fransua/ecolopy.png
           ?branch=master
   :target: http://travis-ci.org/fransua/ecolopy
   :alt: Build Status
