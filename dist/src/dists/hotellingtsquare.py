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
    def cdf(m,p,x):
        if(x>0):
            return sp.betainc(p/2,1/2*(m-p+1),(x*(m-p+1))/(m*(x*(m-p+1)/m+m-p+1)))
        return 0
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
    def kurtosis(m,p):
        if(m-p>7):
            return 3*(m-p-3)*(-(m-5)*p**2+(m+3)*(m-1)*p+4*(m-1)**2)/((m-1)*p*(m-p-7)*(m-p-5))
        return None
    @staticmethod
    def entropy():
        pass
    @staticmethod
    def skewness(m,p):
        if(m-p>5):
            return 2*(m+p-1)*math.sqrt((-2*m+2*p+6)/(p-m*p))/(m-p-5)
        return None
    @staticmethod
    def ppf():
        pass
    @staticmethod
    def mle():
        pass