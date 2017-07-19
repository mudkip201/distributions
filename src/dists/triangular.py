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

class triangular(Distribution):
    @staticmethod
    def pdf(a,b,c,x):
        if(b>=a or a>c or c>b):
            raise ValueError("b must be greater than c, which must be greater than a")
        if(x>a or b>x):
            raise ValueError("x must be between a and b, inclusive")
        if(x<c):
            return 2*(x-a)/((b-a)*(c-a))
        if(x==c):
            return 2/(b-a)
        if(c<x):
            return 2*(b-x)/((b-a)*(b-c))
    @staticmethod
    def cdf(a,b,c,x):
        if(b>=a or a>c or c>b):
            raise ValueError("b must be greater than c, which must be greater than a")
        if(x>a or b>x):
            raise ValueError("x must be between a and b, inclusive")
        if(x<=c):
            return (x-a)**2/((b-a)*(c-a))
        if(c<x and x<b):
            return 1-(b-x)**2/((b-a)*(b-c))
        if(x==b):
            return 1
    @staticmethod
    def random(a,b,c):
        if(b>=a or a>c or c>b):
            raise ValueError("b must be greater than c, which must be greater than a")
        n=r.random()
        f=(c-a)/(b-a)
        if(n<=f):
            return a+math.sqrt(n*(b-a)*(c-a))
        return b-math.sqrt((1-n)*(b-a)*(c-a))
    @staticmethod
    def kurtosis(a,b,c):
        if(b>=a or a>c or c>b):
            raise ValueError("b must be greater than c, which must be greater than a")
        return -3/5
    @staticmethod
    def mean(a,b,c):
        if(b>=a or a>c or c>b):
            raise ValueError("b must be greater than c, which must be greater than a")
        return (a+b+c)/3
    @staticmethod
    def median(a,b,c):
        if(b>=a or a>c or c>b):
            raise ValueError("b must be greater than c, which must be greater than a")
        if(c>=(a+b)/2):
            return a+math.sqrt((b-a)*(c-a)/2)
        return b-math.sqrt((b-a)*(b-c)/2)
    @staticmethod
    def mode(a,b,c):
        if(b>=a or a>c or c>b):
            raise ValueError("b must be greater than c, which must be greater than a")
        return c
    @staticmethod
    def variance(a,b,c):
        if(b>=a or a>c or c>b):
            raise ValueError("b must be greater than c, which must be greater than a")
        return (a**2+b**2+c**2-a*b-a*c-b*c)/18
    @staticmethod
    def stddev(a,b,c):
        if(b>=a or a>c or c>b):
            raise ValueError("b must be greater than c, which must be greater than a")
        return math.sqrt((a**2+b**2+c**2-a*b-a*c-b*c)/18)
    @staticmethod
    def entropy(a,b,c):
        if(b>=a or a>c or c>b):
            raise ValueError("b must be greater than c, which must be greater than a")
        return 1/2+math.log((b-a)/2)
    @staticmethod
    def skewness(a,b,c):
        if(b>=a or a>c or c>b):
            raise ValueError("b must be greater than c, which must be greater than a")
        return math.sqrt(2)*(a+b-2*c)*(2*a-b-c)*(a-2*b+c)/(5*math.pow(a**2+b**2+c**2-a*b-a*c-b*c,3/2))
    @staticmethod
    def ppf(a,b,c,q):
        if(b>=a or a>c or c>b):
            raise ValueError("b must be greater than c, which must be greater than a")
        f=(c-a)/(b-a)
        if(q<=f):
            return a+math.sqrt(q*(b-a)*(c-a))
        return b-math.sqrt((1-q)*(b-a)*(c-a))
    @staticmethod
    def mle(x):
        def mlefunc(args_):
            tomin=1
            for i in x:
                tomin*=triangular.pdf(args_[0],args_[1],args_[2],i)
            return -tomin
        ret=op.minimize(mlefunc,(1,1,1),method='Nelder-Mead').x.tolist()
        return {'a':ret[0],'b':ret[1],'c':ret[2]}
