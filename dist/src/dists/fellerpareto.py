'''
Created on Jul 15, 2017

@author: matthewcowen-green
'''


import dists.Distribution.Distribution as Distribution
import math
import dists.gamma.gamma as gamma

class fellerpareto(Distribution):
    @staticmethod
    def random(mu,sigma,gmma,d1,d2):
        return mu+sigma*math.pow(gamma.random(d1,1)/gamma.random(d2,1),gmma)
    @staticmethod
    def mode(mu,sigma,gmma,d1,d2):
        return mu+sigma*math.pow((d2-gmma)/(d1+gmma),gmma)
