'''
Created on Jul 16, 2017

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

class johnsonsu(Distribution):
    @staticmethod
    def random(delta,gmma,xi,lmbda):
        return lmbda*math.sinh((normal.random(0,1)-gmma)/delta)+xi
    @staticmethod
    def pdf(delta,gmma,xi,lmbda,x):
        return delta/(gmma*math.sqrt(2*math.pi))*1/math.sqrt(1+((x-xi)/lmbda)**2)*math.exp(-1/2*math.pow(gmma+delta*math.asinh((x-xi)/lmbda),2))
    @staticmethod
    def cdf(delta,gmma,xi,lmbda,x):
        return st.norm.cdf(gmma+delta*math.asinh((x-xi)/lmbda))
    @staticmethod
    def kurtosis(delta,gmma,xi,lmbda):
        z=1/delta**2
        return (4*math.exp((2+2*gmma*delta)*z)*(2+math.exp(z))+4*math.exp((2+6*gmma*delta)*z)*(2+math.exp(z))+6*math.exp(4*gmma*z)*(1+2*math.exp(z))+math.exp(2*z)*(-3+math.exp(2*z)*(3+math.exp(z)*(2+math.exp(z))))+math.exp((2+8*gmma*delta)*z)*(-3+math.exp(2*z)*(3+math.exp(z)*(2+math.exp(z)))))/math.pow(math.exp(z)+2*math.exp(2*z)+math.exp((1+4*gmma*delta)*z),2)
    @staticmethod
    def mean(delta,gmma,xi,lmbda):
        return xi-lmbda*math.exp(1/(2*delta**2))*math.sinh(gmma/delta)
    @staticmethod
    def variance(delta,gmma,xi,lmbda):
        return lmbda**2/2*(math.exp(1/delta**2)-1)*(math.exp(1/delta**2)*math.cosh(2*gmma/delta)+1)
    @staticmethod
    def stddev(delta,gmma,xi,lmbda):
        return math.sqrt(lmbda**2/2*(math.exp(1/delta**2)-1)*(math.exp(1/delta**2)*math.cosh(2*gmma/delta)+1))
    @staticmethod
    def skewness(delta,gmma,xi,lmbda):
        i=math.exp(1/(2*delta**2))*math.sqrt(-1+math.exp(1/delta**2))*(3*math.exp(2*gmma/delta)-3*math.exp(4*gmma/delta)+(2+math.exp(1/delta**2)*math.exp(1/delta**2))-(2+math.exp(1/delta**2))*math.exp((1+6*gmma*delta)/delta**2))/math.pow(math.exp(1/delta**2)+2*math.exp(2*gmma/delta)+math.exp((1+4*gmma*delta)/delta**2),3/2)
    @staticmethod
    def mle(x):
        def mlefunc(args_):
            tomin=1
            for i in x:
                tomin*=johnsonsu.pdf(args_[0],args_[1],args_[2],args_[3],i)
            return -tomin
        ret=op.minimize(mlefunc,(1,1,1,1),method='Nelder-Mead').x.tolist()
        return {'delta':ret[0],'gamma':ret[1],'xi':ret[2],'lmbda':ret[3]}