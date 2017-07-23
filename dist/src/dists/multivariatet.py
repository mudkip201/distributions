'''
Created on Jul 23, 2017

@author: matthewcowen-green
'''


import dists.Distribution.Distribution as Distribution
import math
import numpy as np
import scipy.special as sp

class MyClass(Distribution):
    @staticmethod
    def pdf(m,S,nu,x):
        p=S.shape[1]
        ff=sp.gamma((nu+p)/2)/(sp.gamma(nu/2)*math.pow(nu,p/2)*math.pow(math.pi,p/2)*math.sqrt(np.linalg.det(S)))
        ff*math.pow(1+1/nu*np.dot(np.dot((x-m).T,S.I),x-m),-(nu+p)/2)
        pass
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