'''
Created on Jul 16, 2017

@author: matthewcowen-green
'''


import dists.Distribution.Distribution as Distribution
import dists.Distribution as ds
import math
from numpy import random as r
import scipy.optimize as op

class qexp(Distribution): #q-exponential
    @staticmethod
    def random(q1,lmbda):
        return (-(1/(2-q1))*ds.qlog(r.random(),1/(2-q1)))/lmbda
    @staticmethod
    def pdf(q1,lmbda,x):
        return (2-q1)*lmbda*ds.qexp(-lmbda*x,q1)
    @staticmethod
    def cdf(q1,lmbda,x):
        return 1-ds.qexp(-lmbda*x/(1/(2-q1)),1/(2-q1))
    @staticmethod
    def kurtosis(q1,lmbda):
        if(q1<6/5):
            return 6*(-4*q1**3+17*q1**2-20*q1+6)/((q1-2)*(4*q1-5)*(5*q1-6))
        return None
    @staticmethod
    def mean(q1,lmbda):
        if(q1<3/2):
            return 1/(lmbda*(3-2*q1))
        return None
    @staticmethod
    def median(q1,lmbda):
        qp=1/(2-q1)
        return -qp*ds.qlog(1/2,qp)/lmbda
    @staticmethod
    def mode(q1,lmbda):
        return 0
    @staticmethod
    def variance(q1,lmbda):
        if(q1<4/3):
            return (q1-2)/((2*q1-3)**2*(3*q1-4)*lmbda**2)
        return None
    @staticmethod
    def stddev(q1,lmbda):
        if(q1<4/3):
            return math.sqrt((q1-2)/((2*q1-3)**2*(3*q1-4)*lmbda**2))
        return None
    @staticmethod
    def skewness(q1,lmbda):
        if(q1<5/4):
            return 2/(5-4*q1)*math.sqrt((3*q1-4)/(q1-2))
        return None
    @staticmethod
    def ppf(q1,lmbda,q):
        return (-(1/(2-q1))*ds.qlog(q,1/(2-q1)))/lmbda
    @staticmethod
    def mle(x):
        def mlefunc(args_):
            tomin=1
            for i in x:
                tomin*=qexp.pdf(args_[0],args_[1],i)
            return -tomin
        ret=op.minimize(mlefunc,(1,1),method='Nelder-Mead').x.tolist()
        return {'q1':ret[0],'lambda':ret[1]}