'''
Created on Jul 16, 2017

@author: matthewcowen-green
'''


import dists.Distribution.Distribution as Distribution
import dists.Distribution as ds
import math
import scipy.optimize as op

class weibulluniform(Distribution):
    @staticmethod
    def random(a,b,t):
        return t/(1+math.pow(-math.log(1-ds.rg0())/a,-1/b))
    @staticmethod
    def pdf(a,b,t,x):
        return t*a*b/((t-x)**2)*math.pow(x/(t-x),b-1)*math.exp(-a*math.pow(x/(t-x),b))
    @staticmethod
    def cdf(a,b,t,x):
        return 1-math.exp(-a*math.pow(x/(t-x),b))
    @staticmethod
    def median(a,b,t):
        return t/(1+math.pow(-math.log(1/2)/a,-1/b))
    @staticmethod
    def ppf(a,b,t,q):
        return t/(1+math.pow(-math.log(1-q)/a,-1/b))
    @staticmethod
    def mle(x):
        def mlefunc(args_):
            tomin=1
            for i in x:
                tomin*=weibulluniform.pdf(args_[0],args_[1],args_[2],i)
            return -tomin
        ret=op.minimize(mlefunc,(1,1,1),method='Nelder-Mead').x.tolist()
        return {'a':ret[0],'b':ret[1],'t':ret[2]}