#!/usr/bin/python
"""
24 Apr 2012


"""

__author__  = "Francois-Jose Serra"
__email__   = "francois@barrabin.org"
__licence__ = "GPLv3"
__version__ = "0.0"

import os

for root, _, files in os.walk('/home/francisco/Box/untbgen/'):
    print root
    for f in files:
        if not f.endswith('.py'): continue
        lines = open(root + '/' + f).readlines()
        out = open(root + '/' + f, 'w')
        for line in lines:
            if 'import' in line:
                print '- ' + line
                line = line.replace('ecolopy', 'ecolopy_dev')
                print '+ ' + line
            out.write(line)
        out.close()
