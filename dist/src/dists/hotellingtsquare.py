'''
Created on Jul 23, 2017

@author: matthewcowen-green
'''


import dists.Distribution.Distribution as Distribution
import math
import scipy.special as sp

class hotellingtsquare(Distribution):
    @staticmethod
    def pdf(m,p,x):
        return math.pow(x/m,p/2)*math.pow(m/(m+x),(1+m)/2)/(x*sp.beta(p/2,(m-p+1)/2))
    @staticmethod
    def cdf():
        pass
    @staticmethod
    def random():
        pass
    @staticmethod
    def mean(m,p):
        if(m-p-1>0):
            return m*p/(m-p-1)
        return None
    @staticmethod
    def median():
        pass
    @staticmethod
    def mode():
        pass
    @staticmethod
    def variance(m,p):
        if(m-p>3):
            return 2*(m-1)*m**2*p/((m-p-3)*(p-m+1)**2)
        return None
    @staticmethod
    def stddev(m,p):
        if(m-p>3):
            return math.sqrt(2*(m-1)*m**2*p/((m-p-3)*(p-m+1)**2))
        return None
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