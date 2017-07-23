'''
Created on Jul 23, 2017

@author: matthewcowen-green
'''


import dists.Distribution.Distribution as Distribution
import dists.Distribution as ds
import math
import numpy as np
import scipy.linalg as la
import dists.normal.normal as normal

class multivariatenormal(Distribution):
    @staticmethod
    def pdf(m,S,x):
        return math.exp(-1/2*(np.dot(np.dot((x-m).T,S.I),(x-m))))/math.sqrt(2*math.pi*np.linalg.det(S))
    @staticmethod
    def cdf():
        pass
    @staticmethod
    def random(m,S):
        A=np.asmatrix(la.sqrtm(S))
        z=np.zeros(m.shape)
        for i in range(m.shape[1]):
            z[0][i]=normal.random(0,1)
        return m+np.dot(A,z)
    @staticmethod
    def mean(m,S):
        return m
    @staticmethod
    def median():
        pass
    @staticmethod
    def mode(m,S):
        return m
    @staticmethod
    def variance(m,S):
        return S
    @staticmethod
    def stddev():
        pass
    @staticmethod
    def kurtosis():
        pass
    @staticmethod
    def entropy(m,S):
        return 1/2*math.log(np.linalg.det(2*math.pi*ds.euler_gamma*S))
    @staticmethod
    def skewness():
        pass
    @staticmethod
    def ppf():
        pass
    @staticmethod
    def mle(x):
        #x is a set of column vectors x1, x2, x3...
        xx=x.T
        estm=np.average(x,axis=0)
        estS=0
        for i in range(x.dims[0]):
            xi=xx[i].T
            estS+=np.dot((xi-estm).T,(xi-estm))
        estS/=x.dims[0]
        return {'mu':estm,'Sigma':estS}