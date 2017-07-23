'''
Created on Jul 23, 2017

@author: matthewcowen-green
'''


import dists.Distribution.Distribution as Distribution
import math
import scipy.special as sp

class yule(Distribution):
    @staticmethod
    def pdf(a,k):
        return a*sp.beta(1+a,1+k)
    @staticmethod
    def cdf(a,k):
        return 1-a*sp.beta(a,2+math.floor(k))
    @staticmethod
    def random():
        pass
    @staticmethod
    def mean(a):
        if(a>1):
            return 1/(-1+a)
        return None
    @staticmethod
    def median():
        pass
    @staticmethod
    def mode():
        pass
    @staticmethod
    def variance(a):
        if(a>2):
            return -1/(-1+a)**2+a*sp.hyp2f1(2,2,3+a,1)/((-2+a)*(1+a))
        return None
    @staticmethod
    def stddev(a):
        if(a>2):
            return math.sqrt(-1/(-1+a)**2+a*sp.hyp2f1(2,2,3+a,1)/((-2+a)*(1+a)))
        return None
        pass
    @staticmethod
    def kurtosis(a):
        if(a>4):
            return (a-2)*(a*(a*(a+12)-6)+11)/((a-4)*(a-3)*a)
        return None
    @staticmethod
    def entropy():
        pass
    @staticmethod
    def skewness(a):
        if(a>3):
            return math.sqrt(a-2)*(a+1)**2/((a-3)*a)
        return None
    @staticmethod
    def ppf():
        pass
    @staticmethod
    def mle():
        pass