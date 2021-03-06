'''
Created on Jul 19, 2017

@author: matthewcowen-green
'''


import dists.Distribution.Distribution as Distribution
import math
import scipy.special as sp

class genbeta2(Distribution):
    @staticmethod
    def pdf(a,b,p,q,x):
        return a*math.pow(x,a*p-1)/(math.pow(b,a*p)*sp.beta(p,q)*math.pow(1+math.pow(x/b,a),p+q))
    @staticmethod
    def cdf(a,b,p,q,x):
        return math.pow(math.pow(x/b,a)/(1+math.pow(x/b,a)),p)/(p*sp.beta(p,q))*sp.hyp2f1(p,1-q,math.pow(x/b,a))
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