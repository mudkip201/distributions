'''
Created on Jul 16, 2017

@author: matthewcowen-green
'''


import dists.Distribution.Distribution as Distribution
import math
import dists.poisson.poisson as poisson

class hermite(Distribution):
    @staticmethod
    def random(a1,a2):
        return poisson.random(a1)+2*poisson.random(a2)
    @staticmethod
    def kurtosis(a1,a2):
        return (a1+16*a2)/(a1+4*a2)**2
    @staticmethod
    def mean(a1,a2):
        return a1+2*a2
    @staticmethod
    def variance(a1,a2):
        return a1+4*a2
    @staticmethod
    def stddev(a1,a2):
        return math.sqrt(a1+4*a2)
    @staticmethod
    def skewness(a1,a2):
        return (a1+8*a2)/math.pow(a1+4*a2,3/2)
