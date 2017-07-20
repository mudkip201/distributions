'''
Created on Jul 19, 2017

@author: matthewcowen-green
'''


import dists.Distribution.Distribution as Distribution
import math

class discretelogistic(Distribution):
    @staticmethod
    def pdf(m,p,x):
        x=math.floor(x)
        return (1-p)*math.pow(p,x-m)/((1+math.pow(p,x-m))*(1+math.pow(p,x-m+1)))
    @staticmethod
    def cdf(m,p,x):
        return 1/(1+math.pow(p,math.floor(x)-m+1))
    @staticmethod
    def random():
        pass
    @staticmethod
    def mean(m,p):
        return m-0.5
    @staticmethod
    def median():
        pass
    @staticmethod
    def mode():
        pass
    @staticmethod
    def variance(m,p):
        return (math.pi/math.sqrt(3)/math.log(p))**2
    @staticmethod
    def stddev(m,p):
        return math.pi/math.sqrt(3)/math.log(p)
    @staticmethod
    def kurtosis():
        pass
    @staticmethod
    def entropy():
        pass
    @staticmethod
    def skewness():
        pass
    @staticmethod
    def ppf():
        pass
    @staticmethod
    def mle():
        pass