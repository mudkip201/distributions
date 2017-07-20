'''
Created on Jul 15, 2017

@author: matthewcowen-green
'''


import dists.Distribution.Distribution as Distribution
import dists.Distribution as ds
import math
import scipy.optimize as op

class genlogistic(Distribution): #generalized logistic
    @staticmethod
    def random(aa,mu,sigma):
        n=ds.rg0()
        return mu-sigma*math.log((1-math.pow(n,aa))/math.pow(n,aa))
    @staticmethod
    def pdf(aa,mu,sigma,x):
        return aa/math.pow(math.exp((x-mu)/sigma)*(1+math.exp(-(x-mu)/sigma)),aa+1)
    @staticmethod
    def cdf(aa,mu,sigma,x):
        return 1/math.pow(1+math.exp(-(x-mu)/sigma),aa)
    @staticmethod
    def median(aa,mu,sigma):
        return mu-sigma*math.log((1-math.pow(1/2,aa))/math.pow(1/2,aa))
    @staticmethod
    def ppf(aa,mu,sigma,q):
        return mu-sigma*math.log((1-math.pow(q,aa))/math.pow(q,aa))
    @staticmethod
    def mle(x):
        def mlefunc(args_):
            tomin=1
            for i in x:
                tomin*=genlogistic.pdf(args_[0],args_[1],args_[2],i)
            return -tomin
        ret=op.minimize(mlefunc,(1,1,1),method='Nelder-Mead').x.tolist()
        return {'aa':ret[0],'mu':ret[1],'sigma':ret[2]}