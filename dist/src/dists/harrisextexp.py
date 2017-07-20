'''
Created on Jul 16, 2017

@author: matthewcowen-green
'''


import dists.Distribution.Distribution as Distribution
import dists.Distribution as ds
import math
import scipy.optimize as op

class harrisextexp(Distribution): #harris extended exponential
    @staticmethod
    def random(k,l,t):
        return math.pow(l*k,-1)*math.log((1-t)+t*math.pow(1-ds.rg0(),-k))
    @staticmethod
    def pdf(l,k,t,x):
        return l*math.pow(t,1/k)*math.exp(-l*x)/math.pow(1-(1-t)*math.exp(-l*k*x),1+1/k)
    @staticmethod
    def cdf(l,k,t,x):
        return 1-math.pow(t*math.exp(-l*k*x)/(1-(1-t)*math.exp(-l*k*x)),1/k)
    @staticmethod
    def median(k,l,t):
        return math.pow(l*k,-1)*math.log((1-t)+t*math.pow(1/2,-k))
    @staticmethod
    def ppf(l,k,t,q):
        return math.pow(l*k,-1)*math.log((1-t)+t*math.pow(1-q,-k))
    @staticmethod
    def mle(x):
        def mlefunc(args_):
            tomin=1
            for i in x:
                tomin*=harrisextexp.pdf(args_[0],args_[1],args_[2],i)
            return -tomin
        ret=op.minimize(mlefunc,(1,1,1),method='Nelder-Mead').x.tolist()
        return {'l':ret[0],'k':ret[1],'t':ret[2]}