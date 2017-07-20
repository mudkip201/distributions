'''
Created on Jul 16, 2017

@author: matthewcowen-green
'''


import dists.Distribution.Distribution as Distribution
import dists.Distribution as ds
import math
import scipy.optimize as op

class genweibullexp(Distribution): #generalized weibull-exponential
    @staticmethod
    def random(a,c,t,g):
        return -math.log(1-math.pow(1-math.exp(math.pow(-math.log(1-t*ds.rg0())/g,1/a)),1/c))
    @staticmethod
    def pdf(a,c,t,g,x):
        z=1-math.exp(-t*x)
        return c*a/g*t*math.exp(-t*x)*math.pow(z,c-1)/(1-math.pow(z,c))*math.pow(-math.log(1-math.pow(z,c))/g,a-1)*math.pow(math.exp-(-math.log(1-math.pow(1-math.exp(-t*x),c))/g),a)
    @staticmethod
    def cdf(a,c,t,g,x):
        return 1-math.pow(math.exp-(-math.log(1-math.pow(1-math.exp(-t*x),c))/g),a)
    @staticmethod
    def median(a,c,t,g):
        return -math.log(1-math.pow(1-math.exp(math.pow(-math.log(1-t/2)/g,1/a)),1/c))
    @staticmethod
    def ppf(a,c,t,g,q):
        return -math.log(1-math.pow(1-math.exp(math.pow(-math.log(1-t*q)/g,1/a)),1/c))
    @staticmethod
    def mle(x):
        def mlefunc(args_):
            tomin=1
            for i in x:
                tomin*=genweibullexp.pdf(args_[0],args_[1],args_[2],args_[3],i)
            return -tomin
        ret=op.minimize(mlefunc,(1,1,1,1),method='Nelder-Mead').x.tolist()
        return {'a':ret[0],'c':ret[1],'t':ret[2],'g':ret[3]}