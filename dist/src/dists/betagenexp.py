'''
Created on Jul 17, 2017

@author: matthewcowen-green
'''


import dists.Distribution.Distribution as Distribution
import math
import scipy.special as sp
import dists.beta.beta as beta

class betagenexp(Distribution): #Beta generalized exponential
    @staticmethod
    def pdf(a,b,aa,l,x):
        return aa*l/sp.beta(a,b)*math.exp(-l*x)*math.pow(1-math.exp(-l*x),aa*a-1)*math.pow(1-math.pow(1-math.exp(-l*x),aa),b-1)
    @staticmethod
    def cdf(a,b,aa,l,x):
        return sp.betainc(a,b,math.pow(1-math.exp(-l*x),aa))
    @staticmethod
    def random(a,b,aa,l):
        v=beta.random(a,b)
        return -math.log(1-math.pow(v,1/aa))/l
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