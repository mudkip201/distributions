'''
Created on Jul 15, 2017

@author: matthewcowen-green
'''


import dists.Distribution.Distribution as Distribution
import dists.Distribution as ds
import math
import scipy.optimize as op
import scipy.special as sp

class benini(Distribution):
    @staticmethod
    def pdf(aa,bb,sigma,x):
        if(aa<=0 or bb<=0 or sigma<=0 or x<=sigma):
            raise ValueError("aa, bb, and sigma must be greater than 0, and x must be greater than sigma")
        return math.exp(-aa*math.log(x/sigma)-bb*(math.log(x/sigma))**2)*(aa/x+2*bb*math.log(x/sigma)/x)
    @staticmethod
    def cdf(aa,bb,sigma,x):
        if(aa<=0 or bb<=0 or sigma<=0 or x<=sigma):
            raise ValueError("aa, bb, and sigma must be greater than 0, and x must be greater than sigma")
        return 1-math.exp(-aa*math.log(x/sigma)-bb*(math.log(x/sigma))**2)
    @staticmethod
    def median(aa,bb,sigma):
        if(aa<=0 or bb<=0 or sigma<=0):
            raise ValueError("aa, bb, and sigma must be greater than 0, and x must be greater than sigma")
        return sigma*math.exp((-aa+math.sqrt(aa**2+bb*math.log(16)))/(2*bb))
    @staticmethod
    def mean(aa,bb,sigma):
        return math.sqrt(math.pi)*sigma*math.exp((aa-1)**2/(4*bb))*sp.erfc((aa-1)/(2*math.sqrt(bb)))/(2*math.sqrt(bb))+sigma
    @staticmethod
    def variance(aa,bb,sigma):
        ee=sp.erfc((aa-1)/(2*math.sqrt(bb)))
        vv=1/(4*bb)*math.sqrt(math.pi)*sigma**2
        vv*=math.sqrt(math.pi)*(-math.exp((aa-1)**2/(2*bb)))*ee-4*math.sqrt(bb)*math.exp((aa-1)/(4*bb))*ee+4*math.sqrt(bb)*math.exp((aa-2)**2/(4*bb))*sp.erfc((aa-2)/(2*math.sqrt(bb)))
        return vv
    @staticmethod
    def stddev(aa,bb,sigma):
        return math.sqrt(benini.variance(aa,bb,sigma))
    @staticmethod
    def random(aa,bb,sigma):
        n=ds.rg0()
        if(bb==0):
            return math.pow((1/n)/(1/n-1),1/aa)*sigma
        return sigma*math.exp((math.sqrt(aa**2+4*bb*math.log((1/n)/(1/n-1)))+aa)/(2*bb)-aa/bb)
    @staticmethod
    def mle(x):
        def mlefunc(args_):
            tomin=1
            for i in x:
                tomin*=benini.pdf(args_[0],args_[1],args_[2],i)
            return -tomin
        ret=op.minimize(mlefunc,(1,1,1),method='SLSQP',bounds=[(0,None),(0,None),(0,min(x))]).x.tolist()
        return {'aa':ret[0],'bb':ret[1],'sigma':ret[2]}