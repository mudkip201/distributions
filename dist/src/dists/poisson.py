'''
Created on Jul 16, 2017

@author: matthewcowen-green
'''


import dists.Distribution.Distribution as Distribution
import math
from numpy import random as r

class poisson(Distribution):
    @staticmethod
    def pdf(lmbda,k):
        if(lmbda<=0):
            raise ValueError("lambda must be positive")
        if(k%1!=0 or k<0):
            raise ValueError("k must be a non-negative integer")
        return math.pow(lmbda,k)*math.exp(-lmbda)/math.factorial(k)
    @staticmethod
    def random(lmbda):
        if(lmbda<=0):
            raise ValueError("lambda must be positive")
        l=math.exp(-lmbda)
        k=0
        p=1
        while(p>l):
            k+=1
            p*=r.random()
        return int(k)
    @staticmethod
    def kurtosis(lmbda):
        if(lmbda<=0):
            raise ValueError("lambda must be positive")
        return 1/lmbda
    @staticmethod
    def mean(lmbda):
        if(lmbda<=0):
            raise ValueError("lambda must be positive")
        return lmbda
    @staticmethod
    def median(lmbda):
        if(lmbda<=0):
            raise ValueError("lambda must be positive")
        return math.floor(lmbda+1/3-0.02/lmbda)
    @staticmethod
    def mode(lmbda):
        if(lmbda<=0):
            raise ValueError("lambda must be positive")
        return (math.floor(lmbda)-1,math.floor(lmbda))
    @staticmethod
    def variance(lmbda):
        if(lmbda<=0):
            raise ValueError("lambda must be positive")
        return lmbda
    @staticmethod
    def stddev(lmbda):
        if(lmbda<=0):
            raise ValueError("lambda must be positive")
        return math.sqrt(lmbda)
    @staticmethod
    def skewness(lmbda):
        if(lmbda<=0):
            raise ValueError("lambda must be positive")
        return math.sqrt(1/lmbda)
    @staticmethod
    def mle(x):
        return {'lambda':sum(x)/len(x)}