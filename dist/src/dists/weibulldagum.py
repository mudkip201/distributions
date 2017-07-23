'''
Created on Jul 18, 2017

@author: matthewcowen-green
'''


import dists.Distribution.Distribution as Distribution
import dists.Distribution as ds
import math

class weibulldagum(Distribution):
    @staticmethod
    def pdf(b,bb,d,l,x):
        c1=b*l*d*bb*math.pow(x,-d-1)
        c1*=math.pow(1+l*math.pow(x,-d),-b*bb-1)
        c1/=math.pow(1-math.pow(1+l*math.pow(x,-d),-bb),b+1)
        c1*=math.exp(-math.pow(math.pow(1+l*math.pow(x,-d),bb)-1,-b))
        return c1
    @staticmethod
    def cdf(b,bb,d,l,x):
        return 1-math.exp(-math.pow(math.pow(1+l*math.pow(x,-d),bb)-1,-b))
    @staticmethod
    def random(b,bb,d,l):
        u=ds.rg0()
        return math.pow(l,1/d)*math.pow(math.pow(1+math.pow(1-math.log(1-u),-1/b),1/bb)-1,-1/d)
    @staticmethod
    def mean():
        pass
    @staticmethod
    def median(b,bb,d,l):
        return math.pow(l,1/d)*math.pow(math.pow(1+math.pow(1-math.log(1/2),-1/b),1/bb)-1,-1/d)
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
    def ppf(b,bb,d,l,q):
        return math.pow(l,1/d)*math.pow(math.pow(1+math.pow(1-math.log(1-q),-1/b),1/bb)-1,-1/d)
    @staticmethod
    def mle():
        pass