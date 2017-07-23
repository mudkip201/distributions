'''
Created on Jul 23, 2017

@author: matthewcowen-green
'''


import dists.Distribution.Distribution as Distribution
import math
import scipy.special as sp

class benktandergibrat(Distribution):
    @staticmethod
    def pdf(a,b,x):
        return math.exp(-b*math.log(x)**2)*math.pow(x,-2-a)*(-2*b/a+(1+a+2*b*math.log(x))*(1+2*b*math.log(x)/a))
    @staticmethod
    def cdf(a,b,x):
        return 1-math.exp(-b*math.log(x)**2)*math.pow(x,-a-1)*(1+2*b*math.log(x)/a)
    @staticmethod
    def random():
        pass
    @staticmethod
    def mean(a,b):
        return 1+1/a
    @staticmethod
    def median():
        pass
    @staticmethod
    def mode():
        pass
    @staticmethod
    def variance(a,b):
        return -math.sqrt(b)+a*math.exp((a-1)**2/(4*b))*math.sqrt(math.pi)*sp.erfc((a-1)/(2*math.sqrt(b)))/(a**2*math.sqrt(b))
    @staticmethod
    def stddev(a,b):
        return math.sqrt(benktandergibrat.variance(a,b))
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