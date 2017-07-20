'''
Created on Jul 16, 2017

@author: matthewcowen-green
'''


import dists.Distribution.Distribution as Distribution
import dists.Distribution as ds
import math
import scipy.optimize as op

class logisticuniform(Distribution):
    @staticmethod
    def random(l,t):
        return t*(1-math.exp(1-math.pow(1/ds.rg0(),-1/l)))
    @staticmethod
    def pdf(l,t,x):
        return l/(t-x)*math.pow(math.log(t/(t-x)),-l-1)*math.pow(1+math.pow(math.log(t/(t-x)),-l),-2)
    @staticmethod
    def cdf(l,t,x):
        return 1/(1+math.pow(math.log(t/(t-x)),-l))
    @staticmethod
    def median(l,t):
        return t*(1-math.exp(1-math.pow(2,-1/l)))
    @staticmethod
    def ppf(l,t,q):
        return t*(1-math.exp(1-math.pow(1/q,-1/l)))
    @staticmethod
    def mle(x):
        def mlefunc(args_):
            tomin=1
            for i in x:
                tomin*=logisticuniform.pdf(args_[0],args_[1],i)
            return -tomin
        ret=op.minimize(mlefunc,(1,1),method='Nelder-Mead').x.tolist()
        return {'l':ret[0],'t':ret[1]}