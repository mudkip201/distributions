'''
Created on Jul 23, 2017

@author: matthewcowen-green
'''


import dists.Distribution.Distribution as Distribution
import math
import dists.normal.normal as normal

class argus(Distribution): #ARGUS distribution
    @staticmethod
    def pdf(xi,c,x):
        psi=normal.cdf(0,1,xi)-xi*normal.pdf(0,1,xi)-1/2
        return xi**3/(math.sqrt(2*math.pi)*psi)*x/c**2*math.sqrt(1-x**2/c**2)*math.exp(-1/2*xi**2*(1-x**2/c**2))
    @staticmethod
    def cdf(xi,c,x):
        psi1=normal.cdf(0,1,xi*math.sqrt(1-x**2/c**2))-xi*math.sqrt(1-x**2/c**2)*normal.pdf(0,1,xi*math.sqrt(1-x**2/c**2))-1/2
        psi2=normal.cdf(0,1,xi)-xi*normal.pdf(0,1,xi)-1/2
        return 1-psi1/psi2
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
    def mode(xi,c):
        return c/(xi*math.sqrt(2))*math.sqrt((xi**2-2)+math.sqrt(xi**4+4))
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