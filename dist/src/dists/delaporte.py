'''
Created on Jul 23, 2017

@author: matthewcowen-green
'''


import dists.Distribution.Distribution as Distribution
import math
import scipy.special as sp

class delaporte(Distribution):
    @staticmethod
    def pdf(a,b,l,k):
        ff=0
        for i in range(k+1):
            ff+=sp.gamma(a+i)*math.pow(b,i)*math.pow(l,k-i)*math.exp(-l)/(sp.gamma(a)*math.factorial(i)*math.pow(1+b,a+i)*math.factorial(k-i))
        return ff
    @staticmethod
    def cdf(a,b,l,k):
        ff=0
        for i in range(k+1):
            for j in range(i+1):
                ff+=sp.gamma(a+i)*math.pow(b,i)*math.pow(l,j-i)*math.exp(-l)/(sp.gamma(a)*math.factorial(i)*math.pow(1+b,a+i)*math.factorial(j-i))
        return ff
    @staticmethod
    def random():
        pass
    @staticmethod
    def mean(a,b,l):
        return l+a*b
    @staticmethod
    def median():
        pass
    @staticmethod
    def mode(a,b,l): #returns tuple
        z=(a-1)*b+l
        if(z%1==0):
            return (z,z+1)
        return (math.floor(z))
    @staticmethod
    def variance(a,b,l):
        return l+a*b*(1+b)
    @staticmethod
    def stddev(a,b,l):
        return math.sqrt(l+a*b*(1+b))
    @staticmethod
    def kurtosis(a,b,l):
        return (l+3*l**2+a*b*(1+6*l+6*l*b+7*b+12*b**2+6*b**3+3*a*b+6*a*b**2+3*a*b**3))/math.pow(l+a*b*(1+b),2)
    @staticmethod
    def entropy():
        pass
    @staticmethod
    def skewness(a,b,l):
        return (l+a*b*(1+3*b+2*b**2)/math.pow(l+a*b*(1+b),3/2))
    @staticmethod
    def ppf():
        pass
    @staticmethod
    def mle():
        pass