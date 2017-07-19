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

class qweibull(Distribution):
    @staticmethod
    def pdf(a,l,q,x):
        if(1<q and q<2):
            return a*math.pow(l,a)*(2-q)*math.pow(x,a-1)*math.pow(1+(q-1)*math.pow(l*x,a),-1/(q-1))
        elif(q<1):
            return a*math.pow(l,a)*(2-q)*math.pow(x,a-1)*math.pow(1-(1-q)*math.pow(l*x,a),1/(1-q))
    @staticmethod
    def cdf(a,l,q,x):
        if(1<q and q<2):
            return 1-math.pow(1+(q-1)*math.pow(l*x,a),(q-2)/(q-1))
        elif(q<1):
            return 1-math.pow(1-(1-q)*math.pow(l*x,a),(2-q)/(1-q))
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