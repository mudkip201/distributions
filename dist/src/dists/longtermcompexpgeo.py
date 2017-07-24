'''
Created on Jul 24, 2017

@author: matthewcowen-green
'''


import dists.Distribution.Distribution as Distribution
import dists.Distribution as ds
import math

class longtermcompexpgeo(Distribution):
    @staticmethod
    def pdf(l,p,t,x):
        return l*t*(1-p)*math.exp(-l*x)/((1-t)*math.exp(-l*x)+t)**2
    @staticmethod
    def cdf(l,p,t,x):
        return 1-(math.exp(-l*x)*(1-p*t)+p*t)/(math.exp(-l*x)*(1-t)+t)
    @staticmethod
    def random(l,p,t):
        q=ds.rg0()
        return 1/l*math.log((t*(1-p-q)+q)/(t*(1-p-q)))
    @staticmethod
    def mean():
        pass
    @staticmethod
    def median(l,p,t):
        return 1/l*math.log((t*(1/2-p)+1/2)/(t*(1/2-p)))
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
    def ppf(l,p,t,q):
        return 1/l*math.log((t*(1-p-q)+q)/(t*(1-p-q)))
    @staticmethod
    def mle():
        pass