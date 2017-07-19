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

class fisk(Distribution):
    @staticmethod
    def pdf(a,b,c,x):
        if(b<=0 or c<=0):
            raise ValueError("b and c must be bigger than 0")
        if(x<=a):
            raise ValueError("x must be bigger than or equal to a")
        return c/b*math.pow((x-a)/b,c-1)*math.pow(1+math.pow((x-a)/b,c),-2)
    @staticmethod
    def cdf(a,b,c,x):
        if(b<=0 or c<=0):
            raise ValueError("b and c must be bigger than 0")
        if(x<=a):
            raise ValueError("x must be bigger than or equal to a")
        return 1/(1+math.pow((x-a)/b,-c))
    @staticmethod
    def random(a,b,c):
        if(b<=0 or c<=0):
            raise ValueError("b and c must be bigger than 0")
        u=r.random()
        return a+b*math.pow(u/(1-u),1/c)
    @staticmethod
    def kurtosis(a,b,c):
        if(b<=0 or c<=0):
            raise ValueError("b and c must be bigger than 0")
        if(c>4):
            p=math.pi/c
            return -3*(p*b/math.sin(p))**4+12*b*((p*b)**3)/math.sin(p)**2/math.sin(2*p)-12*(p**2)*(b**4)/math.sin(p)/math.sin(3*p)+4**p*(b**4)/math.sin(4*p)
        return None
    @staticmethod
    def mean(a,b,c):
        if(b<=0 or c<=0):
            raise ValueError("b and c must be bigger than 0")
        if(c>1):
            return a+math.pi*b/c/math.sin(math.pi/c)
        return None
    @staticmethod
    def median(a,b,c):
        if(b<=0 or c<=0):
            raise ValueError("b and c must be bigger than 0")
        return a+b
    @staticmethod
    def mode(a,b,c):
        if(b<=0 or c<=0):
            raise ValueError("b and c must be bigger than 0")
        return a+b*math.pow((c-1)/(c+1),1/c)
    @staticmethod
    def variance(a,b,c):
        if(b<=0 or c<=0):
            raise ValueError("b and c must be bigger than 0")
        if(c>2):
            return 2*math.pi*(b**2)/c/math.sin(2*math.pi/c)-(math.pi**2)*(b**2)/(c**2)*(1/math.sin(math.pi/c)**2)
        return None
    @staticmethod
    def stddev(a,b,c):
        if(b<=0 or c<=0):
            raise ValueError("b and c must be bigger than 0")
        if(c>2):
            return math.sqrt(2*math.pi*(b**2)/c/math.sin(2*math.pi/c)-(math.pi**2)*(b**2)/(c**2)*(1/math.sin(math.pi/c)**2))
        return None
    @staticmethod
    def skewness(a,b,c):
        if(b<=0 or c<=0):
            raise ValueError("b and c must be bigger than 0")
        if(c>3):
            return (2*math.pi/c/math.sin(math.pi/c))**3-6*b*(math.pi*b/c)**2/math.sin(math.pi/c)/math.sin(2*math.pi/c)+3*math.pi*(b**3)/c/math.sin(3*math.pi/c)
        return None
    @staticmethod
    def ppf(a,b,c,q):
        if(b<=0 or c<=0):
            raise ValueError("b and c must be bigger than 0")
        return a+b*math.pow(q/(1-q),1/c)
    @staticmethod
    def mle(x):
        def mlefunc(args_):
            tomin=1
            for i in x:
                tomin*=fisk.pdf(args_[0],args_[1],args_[2],i)
            return -tomin
        ret=op.minimize(mlefunc,(min(x)-1,1,1),method='Nelder-Mead').x.tolist()
        return {'a':ret[0],'b':ret[1],'c':ret[2]}