'''
Created on Jul 15, 2017

@author: matthewcowen-green
'''


import dists.Distribution.Distribution as Distribution
import dists.Distribution as ds
import math
import scipy.optimize as op

class expweibullgeo(Distribution): #exponentiated weibull-geometric
    @staticmethod
    def random(a,b,g,t):
        u=ds.rg0()
        return math.pow(-math.log(1-math.pow(u/(1-t*(1-u)),1/a)),1/g)/b
    @staticmethod
    def pdf(a,b,g,t,x):
        (1-t)*a*g*math.pow(b,g)*math.pow(x,g-1)*math.exp(-math.pow(b*x,g))*math.pow(1-math.exp(-math.pow(b*x,g)),a-1)/(1-t*math.pow(1-math.exp(-math.pow(b*x,g)),a))**2
    @staticmethod
    def cdf(a,b,g,t,x):
        return (1-t)*math.pow(1-math.exp(-math.pow(b*x,g)),a)/(1-t*math.pow(1-math.exp(-math.pow(b*x,g)),a))
    @staticmethod
    def median(a,b,g,t):
        return math.pow(-math.log(1-math.pow(1/(2-t),1/a)),1/g)/b
    @staticmethod
    def ppf(a,b,g,t,q):
        return math.pow(-math.log(1-math.pow(q/(1-t*(1-q)),1/a)),1/g)/b
    @staticmethod
    def mle(x):
        def mlefunc(args_):
            tomin=1
            for i in x:
                tomin*=expweibullgeo.pdf(args_[0],args_[1],args_[2],args_[3],i)
            return -tomin
        ret=op.minimize(mlefunc,(1,1,1,1),method='Nelder-Mead').x.tolist()
        return {'a':ret[0],'b':ret[1],'g':ret[2],'t':ret[3]}
