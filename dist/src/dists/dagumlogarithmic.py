'''
Created on Jul 18, 2017

@author: matthewcowen-green
'''


import dists.Distribution.Distribution as Distribution
import math

class dagumlogarithmic(Distribution):
    @staticmethod
    def pdf(b,d,l,t,x):
        C1=math.pow(1-t*math.pow(1+l*math.pow(x,-d),-b),-1)
        return C1*b*l*d*t*math.pow(x,-d-1)*math.pow(1+l*math.pow(x,-d),-b-1)/(-math.log(1-t))
    @staticmethod
    def cdf(b,d,l,t,x):
        C1=-math.log(1-t*math.pow(1+l*math.pow(x,-d),-b))/(-math.log(1-t))
        return C1*(1-t)
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