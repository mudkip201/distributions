'''
Created on Jul 23, 2017

@author: matthewcowen-green
'''


import dists.Distribution.Distribution as Distribution
import math

class marchenkopastur(Distribution):
    @staticmethod
    def pdf():
        pass
    @staticmethod
    def cdf(l,x):
        if(0<l and l<1):
            if(x<=math.pow(math.sqrt(l)-1,2)):
                return 0
            if(x>=math.pow(math.sqrt(l)+1,2)):
                return 1
            if(math.pow(math.sqrt(l)-1,2)<x and x<math.pow(math.sqrt(l)+1,2)):
                return 1/(2*math.pi*l)*(math.pi*l+math.sqrt(-x**2-(l-1)**2+2*x*(1+l))-(1+l)*math.atan((1-x+l)/math.sqrt(4*l-(1-x+l)**2))+(l-1)*math.atan(((l-1)**2-x*(l+1))/((l-1)*math.sqrt(4*l-(l-x+1)**2))))
        if(x<0):
            return 0
        if(x>=(1+math.sqrt(l))**2):
            return 1
        return 1-1/(2*math.pi*l)*(math.pi-math.sqrt(4*l-(l+1-x)**2)+(1+l)*math.atan((l+1-x)/math.sqrt(4*l-(l+1-x)**2))+(l-1)*math.atan((-(l-1)**2+x*(1+l))/((l-1)*math.sqrt(4*l-(l+1-x)**2))))
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
    def mode():
        pass
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