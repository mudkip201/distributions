'''
Created on Jul 16, 2017

@author: matthewcowen-green
'''


import dists.Distribution.Distribution as Distribution
import math
from numpy import random as r
import dists.halflaha.halflaha as halflaha

class laha(Distribution):
    @staticmethod
    def random(mu,sigma2):
        n=r.random()
        if(n<0.5):
            return -halflaha.random(mu,sigma2)
        return halflaha.random(mu,sigma2)
    @staticmethod
    def pdf(mu,sigma2,x):
        return math.pow(sigma2,3/2)*math.sqrt(2)/(math.pi*(sigma2**2+(x-mu)**4))
    @staticmethod
    def mean(mu,sigma2):
        return mu
    @staticmethod
    def variance(mu,sigma2):
        return sigma2
    @staticmethod
    def stddev(mu,sigma2):
        return math.sqrt(sigma2)