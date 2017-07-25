'''
Created on Jul 25, 2017

@author: matthewcowen-green
'''


import dists.Distribution.Distribution as Distribution
import math
import scipy.special as sp
import dists.beta.beta as beta

class beta3(Distribution):
    @staticmethod
    def pdf(a,b,x):
        return math.pow(2,a)*math.pow(x,a-1)*math.pow(1-x,b-1)/(sp.beta(a,b)*math.pow(1+x,a+b))
    @staticmethod
    def cdf():
        pass
    @staticmethod
    def random(a,b):
        u=beta.random(a,b)
        return u/(2-u)
    @staticmethod
    def mean(a,b):
        return a/(math.pow(2,b)*(a+b))*sp.hyp2f1(b,a+b,a+b+1,1/2)
    @staticmethod
    def median():
        pass
    @staticmethod
    def mode():
        pass
    @staticmethod
    def variance(a,b):
        return a/(math.pow(2,b)*(a+b))*((a+1)/(a+b+1)*sp.hyp2f1(b,a+b,a+b+2,1/2)-sp.hyp2f1(b,a+b,a+b+1,1/2)**2)
    @staticmethod
    def stddev(a,b):
        return math.sqrt(beta3.variance(a,b))
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