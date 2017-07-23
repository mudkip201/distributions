'''
Created on Jul 23, 2017

@author: matthewcowen-green
'''


import dists.Distribution.Distribution as Distribution
import math

class poissonconsul(Distribution):
    @staticmethod
    def pdf(l,m,x):
        return math.exp(-x*l-m)*m*math.pow(x*l+m,x-1)/math.factorial(x)
    @staticmethod
    def cdf():
        pass
    @staticmethod
    def random():
        pass
    @staticmethod
    def mean(l,m):
        return m/(1-l)
    @staticmethod
    def median():
        pass
    @staticmethod
    def mode():
        pass
    @staticmethod
    def variance(m,l):
        return m/(1-l)**3
    @staticmethod
    def stddev(m,l):
        return math.sqrt(m/(1-l)**3)
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