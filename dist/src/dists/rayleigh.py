'''
Created on Jul 16, 2017

@author: matthewcowen-green
'''


import dists.Distribution.Distribution as Distribution
import dists.Distribution as ds
import math
import numpy as np

class rayleigh(Distribution):
    @staticmethod
    def pdf(sigma,x):
        if(sigma<=0 or x<0):
            raise ValueError("sigma must be positive and x must be non-negative")
        return x/(sigma**2)*math.exp(-x**2/(2*sigma**2))
    @staticmethod
    def cdf(sigma,x):
        if(sigma<=0 or x<0):
            raise ValueError("sigma must be positive and x must be non-negative")
        return 1-math.exp(-x**2/(2*sigma**2))
    @staticmethod
    def random(sigma):
        if(sigma<=0):
            raise ValueError("sigma must be positive")
        return sigma*math.sqrt(-2*math.log(ds.rg0()))
    @staticmethod
    def kurtosis(sigma):
        if(sigma<=0):
            raise ValueError("sigma must be positive")
        return -(6*math.pi**2-24*math.pi+16)/((4-math.pi)**2)
    @staticmethod
    def mean(sigma):
        if(sigma<=0):
            raise ValueError("sigma must be positive")
        return sigma*math.sqrt(math.pi/2)
    @staticmethod
    def median(sigma):
        if(sigma<=0):
            raise ValueError("sigma must be positive")
        return sigma*math.sqrt(2*math.log(2))
    @staticmethod
    def mode(sigma):
        if(sigma<=0):
            raise ValueError("sigma must be positive")
        return sigma
    @staticmethod
    def variance(sigma):
        if(sigma<=0):
            raise ValueError("sigma must be positive")
        return (4-math.pi)/2*sigma**2
    @staticmethod
    def stddev(sigma):
        if(sigma<=0):
            raise ValueError("sigma must be positive")
        return math.sqrt((4-math.pi)/2)*sigma
    @staticmethod
    def entropy(sigma):
        if(sigma<=0):
            raise ValueError("sigma must be positive")
        return 1+math.log(sigma/math.sqrt(2))+ds.euler_gamma/2
    @staticmethod
    def skewness(sigma):
        if(sigma<=0):
            raise ValueError("sigma must be positive")
        return 2*math.sqrt(math.pi)*(math.pi-3)/math.pow(4-math.pi,3/2)
    @staticmethod
    def ppf(sigma,q):
        if(sigma<=0):
            raise ValueError("sigma must be positive")
        return sigma*math.sqrt(-2*math.log(q))
    @staticmethod
    def mle(x):
        s2=np.sum(np.power(x,2))/(2*len(x))
        n=len(x)
        return math.sqrt(s2)*math.pow(4,n)*math.factorial(n)*math.factorial(n-1)*math.sqrt(n)/(math.factorial(2*n)*math.sqrt(math.pi))