'''
Created on Jul 16, 2017

@author: matthewcowen-green
'''


import dists.Distribution.Distribution as Distribution
import dists.poisson.poisson as poisson
import dists.gamma.gamma as gamma

class polya(Distribution):
    @staticmethod
    def random(a,b):
        return poisson.random(gamma.random(a,b))
