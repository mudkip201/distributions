'''
Created on Jul 16, 2017

@author: matthewcowen-green
'''


import dists.Distribution.Distribution as Distribution
import dists.Distribution as ds
import math
import scipy.optimize as op

class oddgenexpgompertz(Distribution): #odd generalized exponential-gompertz
    @staticmethod
    def random(a,b,c,l):
        return 1/c*math.log(1+c/l*math.log(1-1/a*math.log(1-math.pow(ds.rg0(),1/b))))
    @staticmethod
    def pdf(a,b,c,l,x):
        return a*b*l*math.exp(c*x)*math.exp(l/c*(math.exp(c*x)-1))*math.exp(-a*(math.exp(l/c*(math.exp(c*x)-1))-1))*math.pow(1-math.exp(-a*(math.exp(l/c*(math.exp(c*x)-1))-1)),b-1)
    @staticmethod
    def cdf(a,b,c,l,x):
        return math.pow(1-math.exp(-a*(math.exp(l/c*(math.exp(c*x)-1))-1)),b)
    @staticmethod
    def median(a,b,c,l):
        return 1/c*math.log(1+c/l*math.log(1-1/a*math.log(1-math.pow(1/2,1/b))))
    @staticmethod
    def ppf(a,b,c,l,q):
        return 1/c*math.log(1+c/l*math.log(1-1/a*math.log(1-math.pow(q,1/b))))
    @staticmethod
    def mle(x):
        def mlefunc(args_):
            tomin=1
            for i in x:
                tomin*=oddgenexpgompertz.pdf(args_[0],args_[1],args_[2],args_[3],i)
            return -tomin
        ret=op.minimize(mlefunc,(1,1,1,1),method='Nelder-Mead').x.tolist()
        return {'a':ret[0],'b':ret[1],'c':ret[2],'l':ret[3]}