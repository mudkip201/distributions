'''
Created on Jul 18, 2017

@author: matthewcowen-green
'''


import dists.Distribution.Distribution as Distribution
import dists.Distribution as ds
import math

class burr3poisson(Distribution):
    @staticmethod
    def pdf(b,d,t,x):
        return b*d*t*math.pow(x,-d-1)*math.pow(1+math.pow(x,-d),-b-1)*math.exp(t*math.pow(1+math.pow(x,-d),-b))/(math.exp(t)-1)
    @staticmethod
    def cdf(b,d,t,x):
        return (math.exp(t*math.pow(1+math.pow(x,-d),-b))-1)/(math.exp(t)-1)
    @staticmethod
    def random(b,d,t):
        u=ds.rg0()
        return math.pow((math.pow(math.log(u*(math.exp(t)-1)+1)/t,-1/b)-1),-1/d)
    @staticmethod
    def mean():
        pass
    @staticmethod
    def median(b,d,t):
        return math.pow((math.pow(math.log((math.exp(t)-1)+1)/(2*t),-1/b)-1),-1/d)
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
    def ppf(b,d,t,q):
        return math.pow((math.pow(math.log(q*(math.exp(t)-1)+1)/t,-1/b)-1),-1/d)
    @staticmethod
    def mle():
        pass