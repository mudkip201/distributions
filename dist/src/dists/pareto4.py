'''
Created on Jul 16, 2017

@author: matthewcowen-green
'''


import dists.Distribution.Distribution as Distribution
import math
import dists.fellerpareto.fellerpareto as fellerpareto

class pareto4(Distribution):
    @staticmethod
    def random(mu,sigma,gmma,aa):
        return fellerpareto.random(mu,sigma,gmma,aa,1)
    @staticmethod
    def cdf(mu,sigma,gmma,aa,x):
        if(x>mu):
            return math.pow(1+math.pow((x-mu)/sigma,1/gmma),-aa)
        return 0
