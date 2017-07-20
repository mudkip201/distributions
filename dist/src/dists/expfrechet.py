'''
Created on Jul 15, 2017

@author: matthewcowen-green
'''


import dists.Distribution.Distribution as Distribution
import dists.Distribution as ds
import math
import scipy.optimize as op

class exponentiatedfrechet(Distribution):
    @staticmethod
    def random(a,l,s):
        u=ds.rg0()
        return s/math.pow(-math.log(1-math.pow(1-u,1/a)),1/l)
    @staticmethod
    def pdf(a,l,s,x):
        return a*l*math.pow(s,l)/math.pow(x,-l-1)*math.pow(1-math.exp(-math.pow(s/x,l)),a-1)*math.exp(-math.pow(s/x,l))
    @staticmethod
    def cdf(a,l,s,x):
        return 1-math.pow(1-math.exp(-math.pow(s/x,l)),a)
    @staticmethod
    def median(a,l,s):
        return s/math.pow(-math.log(1-math.pow(1/2,1/a)),1/l)
    @staticmethod
    def ppf(a,l,s,q):
        return s/math.pow(-math.log(1-math.pow(1-q,1/a)),1/l)
    @staticmethod
    def mle(x): #not working
        def mlefunc(args_):
            tomin=1
            for i in x:
                tomin*=exponentiatedfrechet.pdf(args_[0],args_[1],args_[2],i)
            return -tomin
        ret=op.minimize(mlefunc,(1,1,1),method='Nelder-Mead').x.tolist()
        return {'a':ret[0],'l':ret[1],'s':ret[2]}