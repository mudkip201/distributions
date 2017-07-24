'''
Created on Jul 24, 2017

@author: matthewcowen-green
'''


import dists.Distribution.Distribution as Distribution
import math
import dists.uniform.uniform as uniform
import dists.exponential.exponential as exponential
import dists.gamma.gamma as gamma

class quasixgamma(Distribution):
    @staticmethod
    def pdf(a,t,x):
        return t/(1+a)*(a+t**2/2*x**2)*math.exp(-t*x)
    @staticmethod
    def cdf(a,t,x):
        return 1-(1+a+t*x+t**2*x**2/2)/(1+a)*math.exp(-t*x)
    @staticmethod
    def random(a,t):
        u=uniform.random(0,1)
        v=exponential.random(t)
        w=gamma.random(3,t)
        if(u<a/(1+a)):
            return v
        return w
    @staticmethod
    def mean():
        pass
    @staticmethod
    def median():
        pass
    @staticmethod
    def mode(a,t):
        if(0<a and a<=1/2):
            return (1+math.sqrt(1-2*a))/t
        return 0
    @staticmethod
    def variance(a,t):
        return (a**2+8*a+3)/(t**2*(1+a)**2)
    @staticmethod
    def stddev(a,t):
        return math.sqrt((a**2+8*a+3)/(t**2*(1+a)**2))
    @staticmethod
    def kurtosis(a,t):
        return 3*(a**4+88*a**3+310*a**2+288*a+177)/(t**4*(1+a)**4)
    @staticmethod
    def entropy():
        pass
    @staticmethod
    def skewness(a,t):
        return 2*(a**3+15*a**2+9*a+3)/(t**3*(1+a)**3)
    @staticmethod
    def ppf():
        pass
    @staticmethod
    def mle():
        pass