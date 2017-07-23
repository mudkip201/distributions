'''
Created on Jul 23, 2017

@author: matthewcowen-green
'''


import dists.Distribution.Distribution as Distribution
import math
import numpy as np
import scipy.special as sp

class dirichletmultinomial(Distribution):
    @staticmethod
    def pdf(a,n,x):
        ff=math.factorial(n)*sp.gamma(np.sum(a))/sp.gamma(n+np.sum(a))
        for i in range(a.shape[1]):
            ff*=sp.gamma(x[0][i]+a[0][i])/(math.factorial(x[0][i])*sp.gamma(a[0][i]))
        return ff
    @staticmethod
    def cdf():
        pass
    @staticmethod
    def random():
        pass
    @staticmethod
    def mean(a,n):
        return n*a/np.sum(a)
    @staticmethod
    def median():
        pass
    @staticmethod
    def mode():
        pass
    @staticmethod
    def variance(a,n):
        return n*a/np.sum(a)*(1-a/np.sum(a))*((n+np.sum(a))/(1+np.sum(a)))
    @staticmethod
    def covariance(a,n):
        cc=np.dot(a,a.T)
        for i in range(a.shape[1]):
            cc[i][i]=0
        cc/=np.sum(a)**2
        cc*=(n+np.sum(a))/(1+np.sum(a))
        return cc
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