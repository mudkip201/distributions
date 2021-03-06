'''
Created on Jul 23, 2017

@author: matthewcowen-green
'''


import dists.Distribution.Distribution as Distribution
import math

class boreltanner(Distribution):
    @staticmethod
    def pdf(a,n,x):
        return (math.exp(-x*a)*n*math.pow(x,x-n-1)*math.pow(a,x-n))/math.factorial(x-n)
    @staticmethod
    def cdf():
        pass
    @staticmethod
    def random():
        pass
    @staticmethod
    def mean(a,n):
        return n/(1-a)
    @staticmethod
    def median():
        pass
    @staticmethod
    def mode():
        pass
    @staticmethod
    def variance(a,n):
        return n*a/(1-a)**3
    @staticmethod
    def stddev(a,n):
        return math.sqrt(n*a/(1-a)**3)
    @staticmethod
    def kurtosis(a,n):
        return (6*a**2+8*a+1)/(1-a)*a*n+3
    @staticmethod
    def entropy():
        pass
    @staticmethod
    def skewness(a,n):
        return (2*a+1)/math.sqrt((1-a)*a*n)
    @staticmethod
    def ppf():
        pass
    @staticmethod
    def mle():
        pass