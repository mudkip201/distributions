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

class marshallolkinqweibull(Distribution):
    @staticmethod
    def pdf(a,l,p,q,x):
        F=0
        f=0
        if(1<q and q<2):
            F=math.pow(1+(q-1)*math.pow(l*x,a),(q-2)/(q-1))
            f=a*math.pow(l,a)*(2-q)*math.pow(x,a-1)*math.pow(1+(q-1)*math.pow(l*x,a),-1/(q-1))
        elif(q<1):
            F=math.pow(1-(1-q)*math.pow(l*x,a),(2-q)/(1-q))
            f=a*math.pow(l,a)*(2-q)*math.pow(x,a-1)*math.pow(1-(1-q)*math.pow(l*x,a),1/(1-q))
        return p*f/math.pow(1-(1-p*(1-F)),2)
    @staticmethod
    def cdf(a,l,p,q,x):
        F=0
        if(1<q and q<2):
            F=math.pow(1+(q-1)*math.pow(l*x,a),(q-2)/(q-1))
        elif(q<1):
            F=math.pow(1-(1-q)*math.pow(l*x,a),(2-q)/(1-q))
        return p*(1-F)/(1-(1-p)*(1-F))
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