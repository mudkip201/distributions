'''
Created on Jul 18, 2017

@author: matthewcowen-green
'''


import dists.Distribution.Distribution as Distribution
import math

class dagumbinomial(Distribution):
    @staticmethod
    def pdf(b,d,l,m,t,x):
        C1=m*math.pow(1+t*math.pow(1+l*math.pow(x,-d),-b),m-1)
        return C1*b*l*d*t*math.pow(x,-d-1)*math.pow(1+l*math.pow(x,-d),-b-1)/(m*math.pow(1+t,m-1))
    @staticmethod
    def cdf(b,d,l,m,t,x):
        C1=math.pow(1+(t*math.pow(1+l*math.pow(x,-d),-b)),m)-1
        return C1/(math.pow(1+t,m)-1)
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