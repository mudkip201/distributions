'''
Created on Jul 15, 2017

@author: matthewcowen-green
'''


import dists.Distribution.Distribution as Distribution
import math
from numpy import random as r
import scipy.optimize as op

class dagum(Distribution):
    @staticmethod
    def pdf(p,a,b,x):
        if(p<=0 or a<=0 or b<=0 or x<=0):
            raise ValueError("p, a, b, and x must be bigger than 0")
        return a*p/x*math.pow(x/b,a*p)/math.pow(math.pow(x/b,a)+1,p+1)
    @staticmethod
    def cdf(p,a,b,x):
        if(p<=0 or a<=0 or b<=0 or x<=0):
            raise ValueError("p, a, b, and x must be bigger than 0")
        return math.pow(1+math.pow(x/b,-a),-p)
    @staticmethod
    def random(p,a,b):
        if(p<=0 or a<=0 or b<=0):
            raise ValueError("p, a, and b must be bigger than 0")
        return b*math.pow(math.pow(r.random(),-1/p)-1,-1/a)
    @staticmethod
    def mean(p,a,b):
        if a>1:
            return b*math.gamma((a-1)/a)*math.gamma(p+1/a)/math.gamma(p)
        return None
    @staticmethod
    def variance(p,a,b):
        if(a>2):
            return b**2*(math.gamma((a-2)/a)*math.gamma(p)*math.gamma(p+2/a)-math.gamma((a-1)/a)**2*math.gamma(p+1/a)**2)/math.gamma(p)**2
        return None
    @staticmethod
    def stddev(p,a,b):
        if(a>2):
            return math.sqrt(dagum.variance(p,a,b))
        return None
    @staticmethod
    def skewness(p,a,b):
        ga1=math.gamma(1-1/a)
        ga2=math.gamma(1-2/a)
        ga3=math.gamma(1-3/a)
        gp1=math.gamma(p+1/a)
        gp2=math.gamma(p+2/a)
        gp3=math.gamma(p+3/a)
        if(a>3):
            return (2*ga1**3*gp1**3-3*ga2*ga1*math.gamma(p)*gp2*gp1+ga3*math.gamma(p)**2*gp3)/math.pow(ga2*math.gamma(p)*gp2-ga1**2*gp1**2,3/2)
        return None
    @staticmethod
    def kurtosis(p,a,b):
        ga1=math.gamma(1-1/a)
        ga2=math.gamma(1-2/a)
        ga3=math.gamma(1-3/a)
        ga4=math.gamma(1-4/a)
        gp0=math.gamma(p)
        gp1=math.gamma(p+1/a)
        gp2=math.gamma(p+2/a)
        gp3=math.gamma(p+3/a)
        gp4=math.gamma(p+4/a)
        if(a>4):
            return (ga4*gp0**3*gp4+ga1*gp1*(-3*ga1**3*gp1**3+6*ga2*ga1*gp0*gp2*gp1-4*ga3*gp0**2*gp3))/(ga1**2*gp1**2-ga2*gp0*gp2)**2
    @staticmethod
    def median(p,a,b):
        if(p<=0 or a<=0 or b<=0):
            raise ValueError("p, a, and b must be bigger than 0")
        return b*math.pow(-1+math.pow(2,1/p),-1/a)
    @staticmethod
    def ppf(p,a,b,q):
        if(p<=0 or a<=0 or b<=0):
            raise ValueError("p, a, and b must be bigger than 0")
        return b*math.pow(math.pow(q,-1/p)-1,-1/a)
    @staticmethod
    def mle(x):
        def mlefunc(args_):
            tomin=1
            for i in x:
                tomin*=dagum.pdf(args_[0],args_[1],args_[2],i)
            return -tomin
        ret=op.minimize(mlefunc,(3,3,3),method='Nelder-Mead').x.tolist()
        #ret=op.differential_evolution(mlefunc,[(0.01,10),(0.01,10),(0.01,10)]).x.tolist()
        return {'p':ret[0],'a':ret[1],'b':ret[2]}