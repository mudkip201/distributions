'''
Created on Jul 16, 2017

@author: matthewcowen-green
'''


import dists.Distribution.Distribution as Distribution
import math
import dists.normal.normal as normal
import dists.chi2.chi2 as chi2

class noncentralt(Distribution):
    @staticmethod
    def random(mu,nu):
        return (normal.random(0,1)+mu)/math.sqrt(chi2.random(nu)/nu)