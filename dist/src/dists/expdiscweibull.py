'''
Created on Jul 15, 2017

@author: matthewcowen-green
'''


import dists.Distribution.Distribution as Distribution
import dists.Distribution as ds
import math
import scipy.optimize as op

class expdiscweibull(Distribution): #exponentiated discrete weibull
    @staticmethod
    def random(p,a,g):
        return math.pow(math.pow(math.log(1-ds.rg0()),(1/g))/math.log(p),1/a)-1
    @staticmethod
    def pdf(p,a,g,x):
        return math.pow(1-math.pow(p,math.pow(x+1,a)),g) - math.pow(1-math.pow(p,math.pow(x,a)),g)
    @staticmethod
    def cdf(p,a,g,x):
        return math.pow(1-math.pow(p,math.pow(x+1,a)),g)
    @staticmethod
    def median(p,a,g):
        return math.pow(math.pow(math.log(1/2),(1/g))/math.log(p),1/a)-1
    @staticmethod
    def ppf(p,a,g,q):
        return math.pow(math.pow(math.log(1-q),(1/g))/math.log(p),1/a)-1
    @staticmethod
    def mle(x):
        def mlefunc(args_):
            tomin=1
            for i in x:
                tomin*=expdiscweibull.pdf(args_[0],args_[1],args_[2],i)
            return -tomin
        ret=op.minimize(mlefunc,(0.5,1,1),method='Nelder-Mead').x.tolist()
        return {'p':ret[0],'a':ret[1],'g':ret[2]}