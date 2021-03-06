'''
Created on Jul 16, 2017

@author: matthewcowen-green
'''


import dists.Distribution.Distribution as Distribution
import math
import dists.exponential.exponential as exponential
import dists.gamma.gamma as gamma

class k(Distribution):
    @staticmethod
    def random(a,b):
        return math.sqrt(exponential.random(1)*gamma.random(a,b/a))