'''
Created on Jul 16, 2017

@author: matthewcowen-green
'''


import dists.Distribution.Distribution as Distribution
import math
import scipy.optimize as op
import dists.laplace.laplace as laplace

class loglaplace(Distribution):
    @staticmethod
    def random(mu,b):
        return math.exp(laplace.random(mu,b))
    @staticmethod
    def pdf(mu,b,x):
        return 1/(2*b*x)*math.exp(-abs(math.log(x)-mu)/b)
    @staticmethod
    def cdf(mu,b,x):
        if(x>0):
            return (1+(1-math.exp(-abs(math.log(x)-mu)/b))*(abs(math.log(x)-mu)/(math.log(x)-mu)))
    @staticmethod
    def mle(x):
        def mlefunc(args_):
            tomin=1
            for i in x:
                tomin*=loglaplace.pdf(args_[0],args_[1],i)
            return -tomin
        ret=op.minimize(mlefunc,(1,1),method='Nelder-Mead').x.tolist()
        return {'mu':ret[0],'b':ret[1]}