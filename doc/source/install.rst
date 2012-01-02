Download and Install
********************

.. contents::

*Note: for now this section is only relative to installation under Debian-linux*

GNU/Linux
=========

EcoloPy requires python>=2.6 as well as several dependencies:

* **python-scipy** for optimization, and statistical tests.
* **gmpy2** in order to deal with huge numbers necessary to compute K(D,A) in Etienne model

---------------------------------------------------------

Only needed for drawing graph:

* **python-matplotlib**
* **python-numpy**

Install Python libraries:
-------------------------

**Required:**
::

  apt-get install python-scipy

Accessory:

::

  apt-get install python-matplotlib python-scipy


Install GMPY2
-------------

First check your version of python:

::

  python --version

*should be 2.6.x or 2.7.x*

Than according to this install the corresponding dev:

::

  apt-get install python2.x-dev 

and gmp/mpfr libraries:

::

  apt-get install libmpfr-dev libgmp3-dev

and afterward it is necessary to download and install gmpy2 from here:

http://code.google.com/p/gmpy/downloads/list

unzip archive

::
  
  unzip gmpy2-2.xxx.zip
  cd gmpy2-2.xxx
  sudo python setup.py install


Installing EcoloPy
==================

Once done EcoloPy can be downloaded from here:

https://gitorious.org/ecology/ecology/trees/master

Download the package as tar.gz and uncompress it:

::

  tar xzvf ecology-ecology-master.tar.gz

once done, go in the ecology directory and install it:

::

  sudo python setup.py install

finally it is a good thing to test if every thing is working fine.

Go to the test directory and run:

::

  python test_all.py short

