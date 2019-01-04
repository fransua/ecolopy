#!/usr/bin/python
"""
12 Aug 2011

test functionality of abundance class
"""

__author__  = "Francois-Jose Serra"
__email__   = "francois@barrabin.org"
__licence__ = "GPLv3"
__version__ = "0.13"

import unittest

from dataset        import *
from test_community import *
from test_ewens     import *
from test_etienne   import *
from test_lognormal import *

if __name__ == '__main__':
    unittest.main()