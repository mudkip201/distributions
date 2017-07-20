'''
Created on Jul 15, 2017

@author: matthewcowen-green
'''


import dists.Distribution.Distribution as Distribution
import dists.Distribution as ds
import math
import scipy.optimize as op

class extendedgenexp(Distribution): #extended generalized exponential
    @staticmethod
    def random(a,b,l):
        if(b==0):
            return -(1/l)*math.log(1-math.pow(ds.rg0(),1/a))
        return (1/(b*l))*(1-math.pow(1-math.pow(ds.rg0(),1/a),b))
    @staticmethod
    def pdf(a,b,l,x):
        if(b==0):
            return a*l*math.pow(1-math.exp(-l*x),a-1)*math.exp(-l*x)
        return a*l*math.pow(1-math.pow(1-b*l*x,1/b),a-1)*math.pow(1-b*l*x,1/b-1)
    @staticmethod
    def cdf(a,b,l,x):
        if(b==0):
            return math.pow(1-math.exp(-l*x),a)
        return math.pow(1-math.pow(1-b*l*x,1/b),a)
    @staticmethod
    def median(a,b,l):
        if(b==0):
            return -(1/l)*math.log(1-math.pow(1/2,1/a))
        return (1/(b*l))*(1-math.pow(1-math.pow(1/2,1/a),b))
    @staticmethod
    def ppf(a,b,l,q):
        if(b==0):
            return -(1/l)*math.log(1-math.pow(q,1/a))
        return (1/(b*l))*(1-math.pow(1-math.pow(q,1/a),b))
    @staticmethod
    def mle(x):
        def mlefunc(args_):
            tomin=1
            for i in x:
                tomin*=extendedgenexp.pdf(args_[0],args_[1],args_[2],i)
            return -tomin
        ret=op.minimize(mlefunc,(1,1,1),method='Nelder-Mead').x.tolist()
        return {'a':ret[0],'b':ret[1],'l':ret[2]}