'''
Created on Jul 16, 2017

@author: matthewcowen-green
'''


import dists.Distribution.Distribution as Distribution
import dists.Distribution as ds
import math
import scipy.optimize as op

class oddgenexploglog(Distribution): #odd generalized exponential log-logistic
    @staticmethod
    def random(g,l,s,t):
        return s*math.pow(-l*math.log(1-math.pow(ds.rg0(),1/g)),1/t)
    @staticmethod
    def pdf(g,l,s,t,x):
        return g*t/(l*s)*math.pow(x/s,t-1)*math.pow(1-math.exp(-1/l*math.pow(x/s,t)),g-1)*math.exp(-1/l*math.pow(x/s,t))
    @staticmethod
    def cdf(g,l,s,t,x):
        return math.pow(1-math.exp(-1/l*math.pow(x/s,t)),g)
    @staticmethod
    def median(g,l,s,t):
        return s*math.pow(-l*math.log(1-math.pow(1/2,1/g)),1/t)
    @staticmethod
    def mode(g,l,s,t):
        return s*math.pow((1-t)/t*l,1/t)
    @staticmethod
    def ppf(g,l,s,t,q):
        return s*math.pow(-l*math.log(1-math.pow(q,1/g)),1/t)
    @staticmethod
    def mle(x):
        def mlefunc(args_):
            tomin=1
            for i in x:
                tomin*=oddgenexploglog.pdf(args_[0],args_[1],args_[2],args_[3],i)
            return -tomin
        ret=op.minimize(mlefunc,(1,1,1,1),method='Nelder-Mead').x.tolist()
        return {'g':ret[0],'l':ret[1],'s':ret[2],'t':ret[3]}