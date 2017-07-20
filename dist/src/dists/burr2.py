'''
Created on Jul 15, 2017

@author: matthewcowen-green
'''


import dists.Distribution.Distribution as Distribution
import dists.Distribution as ds
import math
import scipy.optimize as op

class burr2(Distribution):
    @staticmethod
    def pdf(r,mu,sigma,x):
        if(r<=0 or sigma<=0):
            raise ValueError("r and sigma must be greater than 0")
        return 1/sigma*r*math.pow(1+math.exp((-x+mu)/sigma),-r-1)/math.exp((x-mu)/sigma)
    @staticmethod
    def cdf(r,mu,sigma,x):
        if(r<=0 or sigma<=0):
            raise ValueError("r and sigma must be greater than 0")
        return 1/math.pow(1+math.exp((x-mu)/sigma),r)
    @staticmethod
    def random(r,mu,sigma): #standard burr ii
        if(r<=0 or sigma<=0):
            raise ValueError("r and sigma must be bigger than 0")
        return mu+sigma*math.log(math.pow(1/ds.rg0(),1/r)-1)
    @staticmethod
    def median(r,mu,sigma):
        return mu+sigma*math.log(math.pow(2,1/r)-1)
    @staticmethod
    def ppf(r,mu,sigma,q):
        if(r<=0 or sigma<=0):
            raise ValueError("r and sigma must be bigger than 0")
        return mu+sigma*math.log(math.pow(1/q,1/r)-1)
    @staticmethod
    def mle(x):
        def mlefunc(args_):
            tomin=1
            for i in x:
                tomin*=burr2.pdf(args_[0],args_[1],args_[2],i)
            return -tomin
        ret=op.minimize(mlefunc,(1,0,1),method='Nelder-Mead').x.tolist()
        return {'r':ret[0],'mu':ret[1],'sigma':ret[2]}