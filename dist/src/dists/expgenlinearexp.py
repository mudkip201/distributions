'''
Created on Jul 15, 2017

@author: matthewcowen-green
'''


import dists.Distribution.Distribution as Distribution
import dists.Distribution as ds
import math
import scipy.optimize as op

class expgenlinearexp(Distribution): #exponentiated generalized linear exponential
    @staticmethod
    def random(a,b,c,d):
        u=ds.rg0()
        if(b!=0):
            return -a/b+1/b*math.sqrt(a**2+2*b*math.pow(-math.log(1-math.pow(u,1/d)),1/c))
        return 1/a*math.pow(-math.log(1-math.pow(u,1/d)),1/c)
    @staticmethod
    def pdf(a,b,c,d,x):
        return c*d*(a+b*x)*math.pow(a*x+b/2*x**2,c-1)*math.pow(1-math.exp(-math.pow(a*x+b/2*x**2,c)),d-1)*math.exp(-math.pow(a*x+b/2*x**2,c))
    @staticmethod
    def cdf(a,b,c,d,x):
        return math.pow(1-math.exp(-math.pow(a*x+b/2*x**2,c)),d)
    @staticmethod
    def median(a,b,c,d):
        if(b!=0):
            return -a/b+1/b*math.sqrt(a**2+2*b*math.pow(-math.log(1-math.pow(1/2,1/d)),1/c))
        return 1/a*math.pow(-math.log(1-math.pow(1/2,1/d)),1/c)
    @staticmethod
    def ppf(a,b,c,d,q):
        if(b!=0):
            return -a/b+1/b*math.sqrt(a**2+2*b*math.pow(-math.log(1-math.pow(q,1/d)),1/c))
        return 1/a*math.pow(-math.log(1-math.pow(q,1/d)),1/c)
    @staticmethod
    def mle(x): #not working
        def mlefunc(args_):
            tomin=1
            for i in x:
                tomin*=expgenlinearexp.pdf(args_[0],args_[1],args_[2],args_[3],i)
            return -tomin
        ret=op.minimize(mlefunc,(1,1,1,1),method='Nelder-Mead').x.tolist()
        return {'a':ret[0],'b':ret[1],'c':ret[2],'d':ret[3]}