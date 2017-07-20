'''
Created on Jul 16, 2017

@author: matthewcowen-green
'''


import dists.Distribution.Distribution as Distribution
import dists.Distribution as ds
import math
import scipy.optimize as op

class kumaraswamyinvweibull(Distribution): #kumaraswamy inverse weibull
    @staticmethod
    def random(a,aa,b,bb):
        return math.pow(-a*aa/math.log(1-math.pow(1-ds.rg0(),1/b)),1/bb)
    @staticmethod
    def pdf(a,aa,b,bb,x):
        return a*b*aa*bb/math.pow(x,bb+1)*math.exp(-aa/math.pow(x,bb))*math.pow(math.exp(-a/math.pow(x,bb)),a-1)*math.pow(1-math.pow(math.exp(-aa/math.pow(x,bb)),a),b-1)
    @staticmethod
    def cdf(a,aa,b,bb,x):
        return 1-math.pow(1-math.exp(-a*aa/math.pow(x,bb)),b)
    @staticmethod
    def median(a,aa,b,bb):
        return math.pow(-a*aa/math.log(1-math.pow(1/2,1/b)),1/bb)
    @staticmethod
    def ppf(a,aa,b,bb,q):
        return math.pow(-a*aa/math.log(1-math.pow(1-q,1/b)),1/bb)
    @staticmethod
    def mle(x):
        def mlefunc(args_):
            tomin=1
            for i in x:
                tomin*=kumaraswamyinvweibull.pdf(args_[0],args_[1],args_[2],args_[3],i)
            return -tomin
        ret=op.minimize(mlefunc,(1,1,1,1),method='Nelder-Mead').x.tolist()
        return {'a':ret[0],'aa':ret[1],'b':ret[2],'bb':ret[3]}