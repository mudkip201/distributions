'''
Created on Jul 16, 2017

@author: matthewcowen-green
'''


import dists.Distribution.Distribution as Distribution
import dists.Distribution as ds
import math
import scipy.optimize as op

class kumaraswamyinvexp(Distribution): #kumaraswamy inverse exponential
    @staticmethod
    def random(a,b,l):
        return -a*l/math.log(1-math.pow(1-ds.rg0(),1/b))
    @staticmethod
    def pdf(a,b,l,x):
        return a*b*l/x**2*math.exp(-a*l/x)*math.pow(1-math.exp(-a*l/x),b-1)
    @staticmethod
    def cdf(a,b,l,x):
        return 1-math.pow(1-math.exp(-a*l/x),b)
    @staticmethod
    def median(a,b,l):
        return -a*l/math.log(1-math.pow(1/2,1/b))
    @staticmethod
    def ppf(a,b,l,q):
        return -a*l/math.log(1-math.pow(1-q,1/b))
    @staticmethod
    def mle(x):
        def mlefunc(args_):
            tomin=1
            for i in x:
                tomin*=kumaraswamyinvexp.pdf(args_[0],args_[1],args_[2],i)
            return -tomin
        ret=op.minimize(mlefunc,(1,1,1),method='Nelder-Mead').x.tolist()
        return {'a':ret[0],'b':ret[1],'l':ret[2]}