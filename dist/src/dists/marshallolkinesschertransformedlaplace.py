'''
Created on Jul 16, 2017

@author: matthewcowen-green
'''


import dists.Distribution.Distribution as Distribution
import math
import scipy.optimize as op
import dists.exponential.exponential as exponential
import dists.normal.normal as normal

class marshallolkinesschertransformedlaplace(Distribution):
    @staticmethod
    def random(b,t):
        w=exponential.random(1/b)
        return 2*t/(1-t**2)*w+math.sqrt(2/(1-t**2))*math.sqrt(w)*normal.random(0,1)
    @staticmethod
    def pdf(b,t,x):
        l=math.sqrt(b*(1-t**2))
        k=l/(t+math.sqrt(l+t**2))
        if(x<0):
            return l*k/(1+k**2)*math.exp(l*x/k)
        return l*k/(1+k**2)*math.exp(-l*k*x)
    @staticmethod
    def cdf(b,t,x):
        l=math.sqrt(b*(1-t**2))
        k=l/(t+math.sqrt(l+t**2))
        if(x<0):
            return k**2/(1+k**2)*math.exp(l*x/k)
        return 1-1/(1+k**2)*math.exp(-l*k*x)
    @staticmethod
    def mean(b,t):
        l=math.sqrt(b*(1-t**2))
        k=l/(t+math.sqrt(l+t**2))
        return (1-k**2)/(l*k)
    @staticmethod
    def variance(b,t):
        l=math.sqrt(b*(1-t**2))
        k=l/(t+math.sqrt(l+t**2))
        return (1+k**4)/(l**2*k**2)
    @staticmethod
    def stddev(b,t):
        l=math.sqrt(b*(1-t**2))
        k=l/(t+math.sqrt(l+t**2))
        return math.sqrt((1+k**4)/(l*k))
    @staticmethod
    def mle(x):
        def mlefunc(args_):
            tomin=1
            for i in x:
                tomin*=marshallolkinesschertransformedlaplace.pdf(args_[0],args_[1],i)
            return -tomin
        ret=op.minimize(mlefunc,(1,1),method='Nelder-Mead').x.tolist()
        return {'b':ret[0],'t':ret[1]}