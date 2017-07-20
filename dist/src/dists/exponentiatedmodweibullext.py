'''
Created on Jul 15, 2017

@author: matthewcowen-green
'''


import dists.Distribution.Distribution as Distribution
import dists.Distribution as ds
import math
import scipy.optimize as op

class exponentiatedmodweibullext(Distribution): #exponentiated modified weibull extension
    @staticmethod
    def random(a,b,g,l): #bad
        return a*math.pow(-math.log(1-1/(a*l)*math.log(1-math.pow(ds.rg0(),1/g))),a/b)
    @staticmethod
    def pdf(a,b,g,l,x):
        return l*b*g*math.pow(x/a,b-1)*math.exp(math.pow(x/a,b)+l*a*(1-math.exp(math.pow(x/a,b))))*math.pow(1-math.exp(l*a*(1-math.exp(math.pow(x/a,b)))),g-1)
    @staticmethod
    def cdf(a,b,g,l,x):
        return math.pow(1-math.exp(l*a*(1-math.exp(math.pow(x/a,b)))),g)
    @staticmethod
    def median(a,b,g,l):
        return a*math.pow(-math.log(1-1/(a*l)*math.log(1-math.pow(1/2,1/g))),a/b)
    @staticmethod
    def ppf(a,b,g,l,q):
        return a*math.pow(-math.log(1-1/(a*l)*math.log(1-math.pow(q,1/g))),a/b)
    @staticmethod
    def mle(x):
        def mlefunc(args_):
            tomin=1
            for i in x:
                tomin*=exponentiatedmodweibullext.pdf(args_[0],args_[1],args_[2],args_[3],i)
            return -tomin
        ret=op.minimize(mlefunc,(1,1,1,1),method='Nelder-Mead').x.tolist()
        return {'a':ret[0],'b':ret[1],'g':ret[2],'l':ret[3]}
