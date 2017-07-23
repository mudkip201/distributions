'''
Created on Jul 15, 2017

@author: matthewcowen-green
'''


import dists.Distribution.Distribution as Distribution
import math
from numpy import random as r


class benford(Distribution):
    @staticmethod
    def pdf(b,x):
        return math.log(1/x+1)/math.log(b)
    @staticmethod
    def cdf(b,x):
        if(1<=x and x<b):
            return math.log(math.floor(x)+1)/math.log(b)
        if(x<1):
            return 0
        return 1
    @staticmethod
    def random(b):
        return 1/(math.pow(b,1-r.random()))
    @staticmethod
    def mean(b):
        return b-math.lgamma(b+1)/math.log(b)
    @staticmethod
    def median(b):
        return math.ceil(math.sqrt(b))-1