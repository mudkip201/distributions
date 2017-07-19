'''
Created on Jul 15, 2017

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
import dists.normal.normal as normal
import dists.exponential.exponential as exponential

class expmodnorm(Distribution): #exponentially modified normal
    @staticmethod
    def random(mu,sigma2,lmbda):
        return normal.random(mu,sigma2)+exponential.random(lmbda)
    @staticmethod
    def pdf(mu,lmbda,sigma2,x):
        return lmbda/2*math.exp(lmbda/2*(2*mu+lmbda*sigma2-2*x))*math.erfc((mu+lmbda*sigma2-x)/(math.sqrt(2*sigma2)))
    @staticmethod
    def cdf(mu,lmbda,sigma2,x):
        u=lmbda*(x-mu)
        v=lmbda*math.sqrt(sigma2)
        return st.norm.cdf(u,0,v)-math.exp(-u+v**2/2+math.log(st.norm.cdf(u,v**2,v)))
    @staticmethod
    def kurtosis(mu,lmbda,sigma2):
        return 3*(1+2/(sigma2*lmbda**2)+3/(lmbda**4*sigma2**2))/math.pow(1+1/(lmbda**2*sigma2),2)-3
    @staticmethod
    def mean(mu,lmbda,sigma2):
        return mu+1/lmbda
    @staticmethod
    def variance(mu,lmbda,sigma2):
        return sigma2+1/lmbda**2
    @staticmethod
    def stddev(mu,lmbda,sigma2):
        return math.sqrt(sigma2+1/lmbda**2)
    @staticmethod
    def skewness(mu,lmbda,sigma2):
        return 2/math.pow(sigma2*lmbda**2,3/2)*math.pow(1+1/(sigma2*lmbda**2),-3/2)
    @staticmethod
    def mle(x):
        def mlefunc(args_):
            tomin=1
            for i in x:
                tomin*=expmodnorm.pdf(args_[0],args_[1],args_[2],i)
            return -tomin
        ret=op.minimize(mlefunc,(1,1,1),method='Nelder-Mead').x.tolist()
        return {'mu':ret[0],'lambda':ret[1],'sigma2':ret[2]}