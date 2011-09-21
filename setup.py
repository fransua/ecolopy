#! /usr/bin/env python

import sys

from setuptools import setup, find_packages

#    ["scipy", "Scipy is only required for the clustering validation functions.", 0],
python_dependencies = [
    ["numpy", "Numpy is required for simple operations.", 0],
    ["scipy", "scipy is required mainly for optimizations of likelihoods.", 0],
    ["gmpy2", "GMPY2 is required for huge numbers generated while getting K(D,A).", 0]
]

TAGS = [
    "Development Status :: 6 - Mature",
    "Environment :: Console",
    "Intended Audience :: Developers",
    "Intended Audience :: Other Audience",
    "Intended Audience :: Science/Research",
    "License :: OSI Approved :: GNU General Public License (GPL)",
    "Natural Language :: English",
    "Operating System :: MacOS",
    "Operating System :: Microsoft :: Windows",
    "Operating System :: POSIX :: Linux",
    "Programming Language :: Python",
    "Topic :: Scientific/Engineering :: Bio-Informatics",
    "Topic :: Software Development :: Libraries :: Python Modules",
    ]

def can_import(mname):
    if mname=="gmpy2":
        try:
            __import__("gmpy2")
        except ImportError:
            return False
        else:
            return True
    elif mname == "scipy":
        try:
            import scipy
        except ImportError:
            return False
        else:
            return True
    elif mname == "numpy":
        try:
            import numpy
        except ImportError:
            return False
        else:
            return True
    else:
        try:
            __import__(mname)
        except ImportError:
            return None
        else:
            return True

def ask(string,valid_values,default=-1,case_sensitive=False):
    """ Asks for a keyborad answer """
    v = None
    if not case_sensitive:
        valid_values = [value.lower() for value in valid_values]
    while v not in valid_values:
        v = raw_input("%s [%s]" % (string,','.join(valid_values)))
        if v == '' and default>=0:
            v = valid_values[default]
        if not case_sensitive:
            v = v.lower()
    return v

print
print "Installing EcoloPy (A python Environment for ecological system exploration)."
print
print "Checking dependencies..."
missing = False
for mname, msg, ex in python_dependencies:
    if not can_import(mname):
        print mname, "cannot be found in your python installation."
        print msg
        missing=True
if missing:
    print "\nHowever, you can still install without such functionality."
    con = ask( "Do you want to continue with the installation anyway?", ["y", "n"])
    if con == "n":
        sys.exit()

# SETUP

ete_version = open("VERSION").readline().strip()
mod_name = ete_version.split("rev")[0]

long_description = open("README").read()
long_description += open("INSTALL").read()
long_description.replace("ecolopy", mod_name)

setup(
    name = mod_name,
    version = ete_version,
    packages = find_packages(),

    requires = [],

    # Project uses reStructuredText, so ensure that the docutils get
    # installed or upgraded on the target machine
    install_requires = [
        ],

    package_data = {
    },
    # metadata for upload to PyPI
    author = "Francois-Jose Serra, and Hernan Dopazo",
    author_email = "serra.francois@gmail.com",
    maintainer = "Francois-Jose Serra",
    maintainer_email = "serra.francois@gmail.com",
    platforms = "OS Independent",
    license = "GPLv3",
    description = "A python Environment for ecological system exploration.",
    long_description = long_description,
    classifiers = TAGS,
    provides = [mod_name],
    keywords = "bioinformatics ecology neutral abundance",
    url = "http://ww.google.com",
    download_url = "http://ww.google.com",
)
