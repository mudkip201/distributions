'''
Created on Jul 15, 2017

@author: matthewcowen-green
'''


import dists.Distribution.Distribution as Distribution
import math
import dists.beta.beta as beta

class betaexp(Distribution): #beta-exponential
    @staticmethod
    def random(a,g):
        return math.log(beta.random(a,g))