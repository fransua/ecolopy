#!/usr/bin/python
"""
21 Sep 2011


"""

__author__  = "Francois-Jose Serra"
__email__   = "francois@barrabin.org"
__licence__ = "GPLv3"
__version__ = "0.13"


class EcologicalModel(object):
    """
    Class representing Ecological models

    :argument name: name of the class, can be either ewens, etienne or lognorm
    :returns: EcologicalModel object
    """
    def __init__(self, community, **kwargs):
        self.community = community
        self._parameters = {}
        self._lnL = float('-inf')

        # old
        #for key in kwargs:
        #    setattr(self, key, kwargs[key])
        #if not 'theta' in kwargs or not 'I' in kwargs:
        #    raise Exception ('Must supply theta and I values\n')


    def __str__(self):
        """
        to print model parameters.
        """
        summary = 'EcologicalModel (object)\n'
        summary += 'Model name                : %s\n' % self.__class__.__name__
        for p in self._parameters:
            summary += '%-24s: %s\n' % (p, self._parameters[p])
        summary += '%-24s: %s\n' % ('lnL', self._lnL)
        return summary

    def random_community(self):
        """
        generic function
        """
        pass
        
    def likelihood(self):
        """
        generic function
        """
        pass
        
    def optimal_params(self):
        """
        generic function
        """
        pass
        
    def get_lnL(self):
        """
        generic function
        """
        return self._lnL
        
    def __set_lnL(self, lnl):
        """
        generic function
        """
        self._lnL = lnl
        
    def doc_lnL(self):
        """
        variable doc
        """
        return "Model likelihood"

    lnL = property(get_lnL, __set_lnL, doc=doc_lnL)
