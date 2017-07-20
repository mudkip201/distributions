'''
Created on Jul 15, 2017

@author: matthewcowen-green
'''


import dists.Distribution.Distribution as Distribution
import dists.Distribution as ds
import math
import scipy.optimize as op

class expgeo(Distribution):
    @staticmethod
    def random(l,t):
        u=ds.rg0()
        return -math.log(t*(1-u)/(t*(1-u)+u))/l
    @staticmethod
    def pdf(l,t,x):
        return l*t*math.exp(-l*x)/math.pow(math.exp(-l*x)*(1-t)+t,2)
    @staticmethod
    def cdf(l,t,x):
        return 1-math.exp(-l*x)/(math.exp(-l*x)*(1-t)+t)
    @staticmethod
    def median(l,t):
        return -math.log(t/(t+1))/l
    @staticmethod
    def mean(l,t):
        return -math.log(t)/(l*(1-t))
    @staticmethod
    def mode(l,t):
        return 1/l*math.log((1-t)/t)
    @staticmethod
    def ppf(l,t,q):
        return -math.log(t*(1-q)/(t*(1-q)+q))/l
    @staticmethod
    def mle(x):
        def mlefunc(args_):
            tomin=1
            for i in x:
                tomin*=expgeo.pdf(args_[0],args_[1],i)
            return -tomin
        ret=op.minimize(mlefunc,(1,1),method='Nelder-Mead').x.tolist()
        return {'l':ret[0],'t':ret[1]}