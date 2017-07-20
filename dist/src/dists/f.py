'''
Created on Jul 15, 2017

@author: matthewcowen-green
'''


import dists.Distribution.Distribution as Distribution
import math
import dists.chi2.chi2 as chi2

class f(Distribution):
    @staticmethod
    def random(d1,d2):
        if(d1<=0 or d2<=0):
            raise ValueError("All inputs must be positive")
        return (chi2.random(d1)/d1)/(chi2.random(d2)/d2)
    @staticmethod
    def kurtosis(d1,d2):
        if(d1<=0 or d2<=0):
            raise ValueError("All inputs must be positive")
        if(d2>8):
            return 12*(d1*(5*d2-22)*(d1+d2-2)+(d2-4)*(d2-2)**2)/(d1*(d2-6)*(d2-8)*(d1+d2-2))
    @staticmethod
    def mean(d1,d2):
        if(d1<=0 or d2<=0):
            raise ValueError("All inputs must be positive")
        if(d2>2):
            return d2/(d2-2)
    @staticmethod
    def mode(d1,d2):
        if(d1<=0 or d2<=0):
            raise ValueError("All inputs must be positive")
        if(d1>2):
            return (d1-2)/d1*d2/(d2+2)
    @staticmethod
    def variance(d1,d2):
        if(d1<=0 or d2<=0):
            raise ValueError("All inputs must be positive")
        if(d2>4):
            return (2*d2**2*(d1+d2-2))/(d1*(d2-2)**2*(d2-4))
    @staticmethod
    def stddev(d1,d2):
        if(d1<=0 or d2<=0):
            raise ValueError("All inputs must be positive")
        if(d2>4):
            return math.sqrt((2*d2**2*(d1+d2-2))/(d1*(d2-2)**2*(d2-4)))
    @staticmethod
    def skewness(d1,d2):
        if(d1<=0 or d2<=0):
            raise ValueError("All inputs must be positive")
        if(d2>6):
            return ((2*d1+d2-2)*math.sqrt(8*(d2-4)))/((d2-6)*math.sqrt(d1*(d1+d2-2)))