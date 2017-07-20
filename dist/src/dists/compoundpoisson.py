'''
Created on Jul 15, 2017

@author: matthewcowen-green
'''


import dists.Distribution.Distribution as Distribution
import dists.poisson.poisson as poisson

class compoundpoisson(Distribution):
    @staticmethod
    def random(lmbda,mu):
        return poisson.random(mu*poisson.random(lmbda))
