'''
Created on Jul 16, 2017

@author: matthewcowen-green
'''


import dists.Distribution.Distribution as Distribution
import math
import scipy.optimize as op
import dists.uniform.uniform as uniform

class logtriangular(Distribution):
    @staticmethod
    def random(a,b,c):
        n=uniform.random(a,c)
        if(math.log(n)<=b):
            return a*math.exp(math.sqrt(n*(math.log(b)-math.log(a))*(math.log(c)-math.log(a))))
        return c*math.exp(-math.sqrt(n*(math.log(c)-math.log(b))*(math.log(c)-math.log(a))))
    @staticmethod
    def pdf(a,b,c,x):
        if(a<x and x<=c):
            return 2*math.log(x/a)/(math.log(b/a)*math.log(c/a))
        if(c<x and x<=b):
            return 2*math.log(b/x)/(math.log(b/a)*math.log(b/c))
    @staticmethod
    def median(a,b,c):
        if(math.log(1/2)<=b):
            return a*math.exp(math.sqrt((math.log(b)-math.log(a))/2*(math.log(c)-math.log(a))))
        return c*math.exp(-math.sqrt((math.log(c)-math.log(b))/2*(math.log(c)-math.log(a))))
    @staticmethod
    def ppf(a,b,c,q):
        if(math.log(q)<=b):
            return a*math.exp(math.sqrt(q*(math.log(b)-math.log(a))*(math.log(c)-math.log(a))))
        return c*math.exp(-math.sqrt(q*(math.log(c)-math.log(b))*(math.log(c)-math.log(a))))
    @staticmethod
    def mle(x):
        def mlefunc(args_):
            tomin=1
            for i in x:
                tomin*=logtriangular.pdf(args_[0],args_[1],args_[2],i)
            return -tomin
        ret=op.minimize(mlefunc,(1,1,1),method='Nelder-Mead').x.tolist()
        return {'a':ret[0],'b':ret[1],'c':ret[2]}