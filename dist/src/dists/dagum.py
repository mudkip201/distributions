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