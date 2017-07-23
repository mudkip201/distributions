'''
Created on Jul 23, 2017

@author: matthewcowen-green
'''


import dists.Distribution.Distribution as Distribution
import math
import scipy.special as sp

class multivarhypergeo(Distribution):
    @staticmethod
    def pdf(c,K,N,n,k):
        ff=1/sp.binom(N,n)
        for i in range(c):
            ff*=sp.binom(K[i][0],k[i][0])
        return ff
    @staticmethod
    def cdf():
        pass
    @staticmethod
    def random():
        pass
    @staticmethod
    def mean(c,K,N,n):
        return n*K/N
    @staticmethod
    def median():
        pass
    @staticmethod
    def mode():
        pass
    @staticmethod
    def variance(c,K,N,n):
        return K/N*(1-K/N)*n*(N-n)/(N-1)
    @staticmethod
    def stddev(c,K,N,n):
        return math.sqrt(multivarhypergeo.variance(c,K,N,n))
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