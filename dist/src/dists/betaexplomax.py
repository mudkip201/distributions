'''
Created on Jul 18, 2017

@author: matthewcowen-green
'''


import dists.Distribution.Distribution as Distribution
import math
import scipy.special as sp

class betaexplomax(Distribution):
    @staticmethod
    def pdf(a,b,l,t,bb,x):
        return bb*t*l/sp.beta(a,b)*math.pow(1+l*x,-t-1)*math.pow(1-math.pow(1+l*x,-t),a*bb-1)*math.pow(1-math.pow(1-math.pow(1+l*x,-t),bb),b-1)
    @staticmethod
    def cdf(a,b,l,t,bb,x):
        return sp.betainc(a,b,math.pow(1-math.pow(1+l*x,-t),bb))
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