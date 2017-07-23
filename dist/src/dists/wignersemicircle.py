'''
Created on Jul 23, 2017

@author: matthewcowen-green
'''


import dists.Distribution.Distribution as Distribution
import math

class wignersemicircle(Distribution):
    @staticmethod
    def pdf(R,x):
        if(x>=-R and x<=R):
            return 2/(math.pi*R**2)*math.sqrt(R**2-x**2)
        return 0
    @staticmethod
    def cdf(R,x):
        if(x<-R):
            return 0
        if(x>R):
            return 1
        return 1/2+x*math.sqrt(R**2-x**2)/(math.pi*R**2)+math.asin(x/R)/math.pi
    @staticmethod
    def random():
        pass
    @staticmethod
    def mean(R):
        return 0
    @staticmethod
    def median(R):
        return 0
    @staticmethod
    def mode(R):
        return 0
    @staticmethod
    def variance(R):
        return R**2/4
    @staticmethod
    def stddev(R):
        return R/2
        pass
    @staticmethod
    def kurtosis(R):
        return 2
    @staticmethod
    def entropy(R):
        return math.log(math.pi*R)-1/2
        pass
    @staticmethod
    def skewness(R):
        return 0
    @staticmethod
    def ppf():
        pass
    @staticmethod
    def mle():
        pass