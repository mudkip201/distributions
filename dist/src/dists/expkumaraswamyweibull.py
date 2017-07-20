'''
Created on Jul 19, 2017

@author: matthewcowen-green
'''


import dists.Distribution.Distribution as Distribution
import math

class expkumaraswamyweibull(Distribution):
    @staticmethod
    def pdf(a,b,c,l,t,x):
        f=t*a*b*c*math.pow(l,c)*math.pow(x,c-1)*math.exp(-math.pow(l*x,c))
        f*=math.pow(1-math.exp(-math.pow(l*x,c)),a-1)*math.pow(1-math.pow(1-math.exp(-math.pow(l*x,c)),a),b-1)*math.pow(1-math.pow(1-math.pow(1-math.exp(-math.pow(l*x,c)),a),b),t-1)
        return f
    @staticmethod
    def cdf():
        pass
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