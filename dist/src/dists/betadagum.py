'''
Created on Jul 18, 2017

@author: matthewcowen-green
'''


import dists.Distribution.Distribution as Distribution
import math
import scipy.special as sp

class betadagum(Distribution):
    @staticmethod
    def pdf(a,b,c,l,d,bb,x):
        return c*bb*l*d*math.pow(x,-d-1)/sp.beta(a,b)*math.pow(1+l*math.pow(x,-d),-a*c*bb-1)*math.pow(1-math.pow(1+l*math.pow(x,-d),-c*bb),b-1)
    @staticmethod
    def cdf(a,b,c,l,d,bb,x):
        return sp.betainc(a,b,math.pow(1+l*math.pow(x,-d),-bb*c))
    @staticmethod
    def random():
        pass
    @staticmethod
    def mean():
        pass
    @staticmethod
    def median():
        pass
    @staticmethod
    def mode():
        pass
    @staticmethod
    def variance():
        pass
    @staticmethod
    def stddev():
        pass
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