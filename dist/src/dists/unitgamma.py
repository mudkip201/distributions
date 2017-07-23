'''
Created on Jul 16, 2017

@author: matthewcowen-green
'''


import dists.Distribution.Distribution as Distribution
import math
import dists.gamma.gamma as gamma

class unitgamma(Distribution):
    @staticmethod
    def random(a,b):
        return math.exp(gamma.random(-1/b,a))