'''
Created on Jul 19, 2017

@author: matthewcowen-green
'''


import dists.Distribution.Distribution as Distribution
import dists.Distribution as ds
import math
import numpy as np
from numpy import random as r
import scipy.stats as st
import scipy.special as sp
import scipy.optimize as op

class modburr3(Distribution):
    @staticmethod
    def pdf(dist,aa,bb,gg,*args,**kwargs):
        G=dist.cdf(*args,**kwargs)
        g=dist.pdf(*args,**kwargs)
        H=G/(1-G)
        return aa*bb*math.pow(H,-bb+1)*math.pow(1+g*math.pow(H,-bb),-aa/gg-1)*g/G**2
    @staticmethod
    def cdf(dist,aa,bb,g,*args,**kwargs):
        return math.pow(1+g*math.pow((dist.cdf(*args,**kwargs))/(1-dist.cdf(*args,**kwargs)),-bb),-aa/g)
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

class kumaraswamy(Distribution):
    @staticmethod
    def pdf(dist,a,b,*args,**kwargs):
        g=dist.pdf(*args,**kwargs)
        G=dist.cdf(*args,**kwargs)
        return a*b*g*math.pow(G,a-1)*math.pow(1-math.pow(G,a),b-1)
    @staticmethod
    def cdf(dist,a,b,*args,**kwargs):
        G=dist.cdf(*args,**kwargs)
        return 1-math.pow(1-math.pow(G,a),b)
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