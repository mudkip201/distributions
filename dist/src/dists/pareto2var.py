'''
Created on Jul 16, 2017

@author: matthewcowen-green
'''


import dists.Distribution.Distribution as Distribution
import dists.Distribution as ds
import math
import scipy.optimize as op

class pareto2var(Distribution): #shifted pareto
    @staticmethod
    def random(a,b,c):
        return a*(1/(math.pow(ds.rg0(),b))-1)+c
    @staticmethod
    def pdf(a,b,c,x):
        return b/a*math.pow(a/(x+a-c),b+1)
    @staticmethod
    def cdf(a,b,c,x):
        return 1-math.pow(a/(x+a-c),b)
    @staticmethod
    def kurtosis(a,b,c):
        return 3*a**4*b*(3*b**2+b+2)/((b-4)*(b-3)*(b-2)*(b-1)**4)
    @staticmethod
    def mean(a,b,c):
        return a/(b-1)+c
    @staticmethod
    def median(a,b,c):
        return a*(math.pow(2,1/b)-1)+c
    @staticmethod
    def mode(a,b,c):
        return a
    @staticmethod
    def variance(a,b,c):
        return a**2*b/((b-2)*(b-1)**2)
    @staticmethod
    def stddev(a,b,c):
        return math.sqrt(a**2*b/((b-2)*(b-1)**2))
    @staticmethod
    def skewness(a,b,c):
        return 2*a**3*b*(b+1)/((b-3)*(b-2)*(b-1)**3)
    @staticmethod
    def ppf(a,b,c,q):
        return a*(1/(math.pow(q,b))-1)+c
    @staticmethod
    def mle(x):
        def mlefunc(args_):
            tomin=1
            for i in x:
                tomin*=pareto2var.pdf(args_[0],args_[1],args_[2],i)
            return -tomin
        ret=op.minimize(mlefunc,(1,1,1),method='Nelder-Mead').x.tolist()
        return {'a':ret[0],'b':ret[1],'c':ret[2]}