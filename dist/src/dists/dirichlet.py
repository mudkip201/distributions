'''
Created on Jul 23, 2017

@author: matthewcowen-green
'''


import dists.Distribution.Distribution as Distribution
import math
import numpy as np
import scipy.special as sp
import dists.gamma.gamma as gamma

class dirichlet(Distribution):
    @staticmethod
    def pdf(K,a,x):
        ff=sp.gamma(np.sum(a))
        for i in range(K):
            ff*=math.pow(x[i][0],a[i][0])/sp.gamma(a[i][0])
        return ff
    @staticmethod
    def cdf():
        pass
    @staticmethod
    def random(a,K):
        yy=np.zeros(a.shape)
        for i in range(K):
            yy[i][0]=gamma.random(a[i][0],1)
        yy/=np.sum(yy)
        return yy
    @staticmethod
    def mean(a,K):
        return a/np.sum(a)
    @staticmethod
    def median():
        pass
    @staticmethod
    def mode():
        pass
    @staticmethod
    def variance(a,K):
        a0=np.sum(a)
        return a*(a0-a)/(a0**2*(a0+1))
    @staticmethod
    def covariance(a,K):
        a0=np.sum(a)
        cc=-np.dot(a,a.T)
        for i in range(K):
            cc[i][i]=0
        cc/=a0**2*(a0+1)
        pass
    @staticmethod
    def stddev():
        pass
    @staticmethod
    def kurtosis():
        pass
    @staticmethod
    def entropy(a,K):
        a0=np.sum(a)
        ee=sp.gamma(np.sum(a))
        for i in range(K):
            ee/=sp.gamma(a[i][0])
        ee=math.log(ee)
        ee-=(K-a0)*sp.digamma(a0)
        for i in range(K):
            ee-=(a[i][0]-1)*sp.digamma(a[i][0])
        return ee
    @staticmethod
    def skewness():
        pass
    @staticmethod
    def ppf():
        pass
    @staticmethod
    def mle():
        pass