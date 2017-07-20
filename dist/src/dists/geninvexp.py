'''
Created on Jul 15, 2017

@author: matthewcowen-green
'''


import dists.Distribution.Distribution as Distribution
import dists.Distribution as ds
import math
import scipy.optimize as op

class geninvexp(Distribution):
    @staticmethod
    def random(a,l):
        return -l/math.log(1-math.pow(1-ds.rg0(),1/a))
    @staticmethod
    def pdf(a,l,x):
        return (a*l/x**2)*math.exp(-l/x)*math.pow(1-math.exp(-l/x),a-1)
    @staticmethod
    def cdf(a,l,x):
        return 1-math.pow(1-math.exp(-l/x),a)
    @staticmethod
    def median(a,l):
        return -l/math.log(1-math.pow(1/2,1/a))
    @staticmethod
    def ppf(a,l,q):
        return -l/math.log(1-math.pow(1-q,1/a))
    @staticmethod
    def mle(x):
        def mlefunc(args_):
            tomin=1
            for i in x:
                tomin*=geninvexp.pdf(args_[0],args_[1],i)
            return -tomin
        ret=op.minimize(mlefunc,(1,1),method='Nelder-Mead').x.tolist()
        return {'a':ret[0],'l':ret[1]}